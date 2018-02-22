/*
 * VisualNRTracker.cc
 * Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
 * Date 13/01/2013
 * 
 * This file implements VisualNRTracker class, which provides tracking 
 * facilities for visual data, instead of synthetic, as NRTracker does.
 * 
 */

#include "EstimationUtils.h"
#include "VisualNRTrackerV1.h"
#include "ATANCamera.h"
#include "DBases.h"
#include "OpenGL.h"
#include "NRMapPoint.h"
#include "NRTrackerData.h"
#include "Interpolation.h"
#include "DeformablePatch.h"

#include <cvd/vision.h>

#include <gvars3/instances.h>
#include <gvars3/GStringUtil.h>

#include <omp.h>

#include <sys/time.h>

using namespace std;
using namespace CVD;
using namespace GVars3;

/**
 * Class constructor, as in Tracker class, it needs the image size, a camera, 
 * a map and a deformation base
 * @param irVideoSize
 * @param c
 * @param m
 * @param b
 */
VisualNRTrackerV1::VisualNRTrackerV1(ImageRef irVideoSize, const ATANCamera &c, 
                                       Map &m, DBases & base):
                                       Tracker(irVideoSize, c, m),mBase(base),
                                       mpDefCoefs(NULL), mpDefCoefs_rigid(NULL),
                                       mpDefCoefs_default(NULL)
{
  string basedisk = GV3::get<string>("bases","");
  unsigned int nbases = getNBasis();

  Matrix<3> R;
  R[0] = makeVector(GV3::get<double>("r00",1),GV3::get<double>("r01",0),GV3::get<double>("r02",0));
  R[1] = makeVector(GV3::get<double>("r10",0),GV3::get<double>("r11",1),GV3::get<double>("r12",0));
  R[2] = makeVector(GV3::get<double>("r20",0),GV3::get<double>("r21",0),GV3::get<double>("r22",1));

  Vector<3> T = makeVector(GV3::get<double>("t0",0),GV3::get<double>("t1",0),GV3::get<double>("t2",0));

  mse3CamFromWorld.get_rotation() = R;
  mse3CamFromWorld.get_translation() = T;
  mse3CamFromWorldPrev = mse3CamFromWorld;

  if(!basedisk.empty())
  {
    bool normalizedModel = GV3::get<bool>("normalizedModel",false);
    cout << "Normalizing the bases? " << normalizedModel << endl;

    mBase.normalizeBases(normalizedModel);

    mpDefCoefs = new double[nbases];
    mpDefCoefs_rigid = new double[nbases];
    //allocating memory for coefficient boostrap
    mpDefCoefs_default = new double[nbases];

    //commonsense initialization of the coeficients.
    //Assume the first shape to be the most significant
    //coeficient set for NR experiment 1
    Vector<Dynamic,double,Reference> vDefCoefs(mpDefCoefs,nbases);
    Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpDefCoefs_rigid,nbases);
    Vector<Dynamic,double,Reference> vDefCoefs_default(mpDefCoefs_default,nbases);
    vDefCoefs = Zeros;
    vDefCoefs_rigid = Zeros;
    vDefCoefs_default = Zeros;

    string coeffsdisk = GV3::get<string>("coeffs","");
    if(!coeffsdisk.empty())
    {
        areCoeffsLoadedFromDisk=true;
        cout << "Coeffs successfully loaded from file" << endl;
    }

    mStateEstimator.setCoefficients(mpDefCoefs, mpDefCoefs_rigid, nbases);
    mStateEstimator.setBases(&mBase);
    //mStateEstimator.rigidShapeCoeffs();

    string interp = GV3::get<string>("interpolation","");
    mpInterpolation=NULL;
    if(!interp.empty())
    {
      if(interp == "NN")
        mpInterpolation = new Interpolation_NN(mCamera, mBase, mpDefCoefs,
                                               mpDefCoefs_rigid, mse3CamFromWorld);
      else if(interp == "LN")
        mpInterpolation = new Interpolation_LN(mCamera, mBase, mpDefCoefs,
                                               mpDefCoefs_rigid, mse3CamFromWorld);
      else
        mpInterpolation = NULL;
    }

    mpInterpolation->compute3DShape();
    mpInterpolation->compute2DShape();

    cout << "Data filled up successfully" << endl;
  }
}

// class destructor
VisualNRTrackerV1::~VisualNRTrackerV1()
{
  if(getNBasis())
  {
    if(mpDefCoefs) delete [] mpDefCoefs;
    if(mpDefCoefs_rigid) delete [] mpDefCoefs_rigid;
    if(mpDefCoefs_default) delete [] mpDefCoefs_default;
    if(mpInterpolation) delete mpInterpolation;
  }
}

// Resets the internal class state
void VisualNRTrackerV1::Reset()
{
  Tracker::Reset();
  //leave the rigid coefficients as they were
  for(unsigned int i=0; i<getNBasis();i++)
  {
    mpDefCoefs[i] = mpDefCoefs_default[i];
  }
}

/**
 * Debug function. Introduce the reference coefficients as one of the imputs
 * 
 * @param imFrame
 * @param coeff
 * @param bDraw
 */
void VisualNRTrackerV1::TrackFrame(CVD::Image<CVD::byte> &imFrame, vector<double> &coeff, bool bDraw)
{
    for(unsigned int i = 0; i < coeff.size(); i++)
    {
        mpDefCoefs[i] = coeff[i] - mpDefCoefs_rigid[i];
    }
    Vector<Dynamic,double,Reference> vDefCoefs(mpDefCoefs, mBase.getNBasis());
    Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpDefCoefs_rigid, mBase.getNBasis());
    TrackFrame(imFrame,bDraw);
}

/**
 * This function tracks the frame, even if there is not a map already
 * available
 * 
 * @param imFrame
 * @param bDraw
 */
void VisualNRTrackerV1::TrackFrame(Image<byte> &imFrame, bool bDraw)
{
  mbDraw = bDraw;
  mMessageForUser.str(""); // Wipe the user message clean

  // Take the input video image, and convert it into the tracker's keyframe struct
  // This does things like generate the image pyramid and find FAST corners
  mCurrentKF.mMeasurements.clear();

  mCurrentKF.MakeKeyFrame_Lite(imFrame);

  if(mnFrame==0)
  {
    alignMapFromDisk();
  }
  // Update the small images for the rotation estimator
  static gvar3<double> gvdSBIBlur("Tracker.RotationEstimatorBlur", 0.75, SILENT);
  static gvar3<int> gvnUseSBI("Tracker.UseRotationEstimator", 1, SILENT);
  static gvar3<int> gvnMinKFIntervalPose("Tracker.MinKFIntervalPose",20,SILENT);
  static gvar3<int> gvnMinKFIntervalDef("Tracker.MinKFIntervalDef",5,SILENT);
  static gvar3<double> gvnreprThreshold("Tracker.ReprojectionThreshold",0.03,SILENT);

  mbUseSBIInit = *gvnUseSBI;
  if(!mpSBIThisFrame)
  {
    mpSBIThisFrame = new SmallBlurryImage(mCurrentKF, *gvdSBIBlur);
    mpSBILastFrame = new SmallBlurryImage(mCurrentKF, *gvdSBIBlur);
  }
  else
  {
    delete  mpSBILastFrame;
    mpSBILastFrame = mpSBIThisFrame;
    mpSBIThisFrame = new SmallBlurryImage(mCurrentKF, *gvdSBIBlur);
  }

  // From now on we only use the keyframe struct!
  mnFrame++;
  cout << mnFrame << endl;

  if(mnFrame == iniframe1 || mnFrame == iniframe2)
    mbUserPressedSpacebar = true;

  if(mbDraw)
  {
    glDrawPixels(mCurrentKF.aLevels[0].im);
    if(GV2.GetInt("Tracker.DrawFASTCorners",0, SILENT))
    {
      glColor3f(1,0,1);  glPointSize(1); glBegin(GL_POINTS);
      for(unsigned int i=0; i<mCurrentKF.aLevels[0].vCorners.size(); i++)
        glVertex(mCurrentKF.aLevels[0].vCorners[i]);
      glEnd();
    }
  }

  if(mMap.IsGood())
  {
    if(mnLostFrames < NUM_LOST_FRAMES)  // .. but only if we're not lost!
    {
      if(mbUseSBIInit)
        CalcSBIRotation();
      ApplyMotionModel();       //
      TrackMap();               //  These three lines do the main tracking work.
      UpdateMotionModel();      //

      AssessTrackingQuality();  //  Check if we're lost or if tracking is poor.

      { // Provide some feedback for the user:
        mMessageForUser << "Tracking Map, quality ";
        if(mTrackingQuality == GOOD)
        {
          mMessageForUser << "good.";
        }
        if(mTrackingQuality == DODGY)
        {
          mMessageForUser << "poor.";
        }
        if(mTrackingQuality == BAD)
        {
          mMessageForUser << "bad.";
        }
        mMessageForUser << " Found:";
        for(int i=0; i<mCurrentKF.MAXLEVELS; i++)
          mMessageForUser << " " << manMeasFound[i] << "/" << manMeasAttempted[i];
          mMessageForUser << " Map: " << mMap.vpPoints.size() << "P, " << mMap.vpKeyFrames.size() << "KF";
      }
    }
    else  // what if there is a map, but tracking has been lost?
    {
      mMessageForUser << "** Attempting recovery **.";
      //first clean up the trails
      mMapTrail.clear();

      if(AttemptRecovery())
      {
        TrackMap();
        AssessTrackingQuality();
      }
    }
    if(mbDraw)
      RenderGrid();
  }
  else // If there is no map, try to make one.
    TrackForInitialMap();

  // GUI interface
  while(!mvQueuedCommands.empty())
  {
    GUICommandHandler(mvQueuedCommands.begin()->sCommand, mvQueuedCommands.begin()->sParams);
    mvQueuedCommands.erase(mvQueuedCommands.begin());
  }
}

/**
 * In case a map is not present yet and it is wanted to be initialized from
 * visual data, this function performs an initial reconstruction based on a
 * virtual stereo pair. It is not recommended to be used in this context as
 * the deformation bases can be not correctly matched with the detected points.
 */
void VisualNRTrackerV1::TrackForInitialMap()
{
  // MiniPatch tracking threshhold.
  static gvar3<int> gvnMaxSSD("Tracker.MiniPatchMaxSSD", 100000, SILENT);
  MiniPatch::mnMaxSSD = *gvnMaxSSD;

  // What stage of initial tracking are we at?
  if(mnInitialStage == TRAIL_TRACKING_NOT_STARTED)
  {
    if(mbUserPressedSpacebar) // First spacebar = this is the first keyframe
    {
      mbUserPressedSpacebar = false;
      TrailTracking_Start();
      mnInitialStage = TRAIL_TRACKING_STARTED;
    }
    else
      mMessageForUser << "Point camera at planar scene and press spacebar to start tracking for initial map." << endl;
    return;
  }

  if(mnInitialStage == TRAIL_TRACKING_STARTED)
  {
    int nGoodTrails = TrailTracking_Advance();  // This call actually tracks the trails
    if(nGoodTrails < 10) // if most trails have been wiped out, no point continuing.
    {
      Reset();
      return;
    }

    // If the user pressed spacebar here, use trails to run stereo and make the intial map..
    if(mbUserPressedSpacebar)
    {
      mbUserPressedSpacebar = false;
      vector<pair<ImageRef, ImageRef> > vMatches; // This is the format the mapmaker wants for the stereo pairs
      for(list<Trail>::iterator i = mlTrails.begin(); i!=mlTrails.end(); i++)
      {
        vMatches.push_back(pair<ImageRef, ImageRef>(i->irInitialPos,
                                                    i->irCurrentPos));
      }
      mnInitialStage = TRAIL_TRACKING_COMPLETE;
      mnLastKeyFrameDropped = mnFrame;
      mlTrails.clear();
    }
    else
      mMessageForUser << "Translate the camera slowly sideways, and press spacebar again to perform stereo init." << endl;
  }
}

/**
 * This function is the one that actually tracks the map when it already exists
 */
void VisualNRTrackerV1::TrackMap()
{
  // Some accounting which will be used for tracking quality assessment:
  for(int i=0; i<mCurrentKF.MAXLEVELS; i++)
    manMeasAttempted[i] = manMeasFound[i] = 0;

  // The Potentially-Visible-Set (PVS) is split into pyramid levels.
  vector<NRTrackerData*> avPVS[mCurrentKF.MAXLEVELS];
  for(int i=0; i<mCurrentKF.MAXLEVELS; i++)
    avPVS[i].reserve(500);

  // For all points in the map..
  for(unsigned int i=0; i<mMap.vpPoints.size(); i++)
  {
    NRMapPoint &p= *(mMap.vpPoints[i]);
    // Ensure that this map point has an associated NRTrackerData struct.
    if(!p.pTData)
    {
        p.pTData = new NRTrackerData(&p,&mBase,i,mpDefCoefs);
        p.pTData->v2Found = Zeros;
        p.pTData->updateDef3D = false;
        p.pTData->bOutlier = false;
        p.pTData->v3Rigid = p.v3WorldPos;
    }
    NRTrackerData &TData = *p.pTData;

    // Project according to current view, and if it's not in the image, skip.
    TData.Project(mse3CamFromWorld*mBase.mse3RelativePose, mCamera);
    if(!TData.bInImage)
	  continue;

    // Calculate camera projection derivatives of this point.
    TData.GetDerivsUnsafe(mCamera);

    // And check what the PatchFinder (included in NRTrackerData) makes of the mappoint in this view..
    TData.nSearchLevel = TData.Finder.CalcSearchLevelAndWarpMatrix(TData.Point, mse3CamFromWorld, TData.m2CamDerivs);
    if(TData.nSearchLevel == -1)
      continue;   // a negative search pyramid level indicates an inappropriate warp for this view, so skip.

    // Otherwise, this point is suitable to be searched in the current image! Add to the PVS.
    TData.bSearched = false;
    TData.bFound = false;
    TData.v2Error_CovScaled = TData.dSqrtInvNoise*(TData.v2Found-TData.v2Image);
    avPVS[TData.nSearchLevel].push_back(&TData);
  }

  // Next: A large degree of faffing about and deciding which points are going to be measured!
  // First, randomly shuffle the individual levels of the PVS.
  // comment this for debugging
  for(int i=0; i<mCurrentKF.MAXLEVELS; i++)
    random_shuffle(avPVS[i].begin(), avPVS[i].end());

  // The next two data structs contain the list of points which will next
  // be searched for in the image, and then used in pose update.
  vector<NRTrackerData*> vNextToSearch;
  vIterationSet.clear();

  // Tunable parameters to do with the coarse tracking stage:
  static gvar3<unsigned int> gvnCoarseMin("Tracker.CoarseMin", 20, SILENT); // Min number of large-scale features for coarse stage
  static gvar3<unsigned int> gvnCoarseMax("Tracker.CoarseMax", 60, SILENT); // Max number of large-scale features for coarse stage
  static gvar3<unsigned int> gvnCoarseRange("Tracker.CoarseRange", 30, SILENT); // Pixel search radius for coarse features
  static gvar3<int> gvnCoarseSubPixIts("Tracker.CoarseSubPixIts", 8, SILENT); // Max sub-pixel iterations for coarse features
  static gvar3<int> gvnCoarseDisabled("Tracker.DisableCoarse", 0, SILENT); // Set this to 1 to disable coarse stage (except after recovery)
  static gvar3<double> gvdCoarseMinVel("Tracker.CoarseMinVelocity", 0.006, SILENT);  // Speed above which coarse stage is used.
  static gvar3<double> gvdErrTh("NRTracker.errorThreshold", 0.03, SILENT);

  double err_th = *gvdErrTh;

  unsigned int nCoarseMax = *gvnCoarseMax;
  unsigned int nCoarseRange = *gvnCoarseRange;

  mbDidCoarse = false;

  // Set of heuristics to check if we should do a coarse tracking stage.
  bool bTryCoarse = true;
  if(*gvnCoarseDisabled ||
     mdMSDScaledVelocityMagnitude < *gvdCoarseMinVel  ||
     nCoarseMax == 0)
    bTryCoarse = false;
  if(mbJustRecoveredSoUseCoarse)
  {
    bTryCoarse = true;
    nCoarseMax *=2;
    nCoarseRange *=2;
  }

  // If we do want to do a coarse stage, also check that there's enough high-level
  // PV map points. We use the lowest-res two pyramid levels (MAXLEVELS-1 and MAXLEVELS-2),
  // with preference to MAXLEVELS-1.
  if(bTryCoarse && avPVS[mCurrentKF.MAXLEVELS-1].size() + avPVS[mCurrentKF.MAXLEVELS-2].size() > *gvnCoarseMin )
  {
    // Now, fill the vNextToSearch struct with an appropriate number of
    // NRTrackerDatas corresponding to coarse map points! This depends on how many
    // there are in different pyramid levels compared to CoarseMin and CoarseMax.

    if(avPVS[mCurrentKF.MAXLEVELS-1].size() <= nCoarseMax)
    { // Fewer than CoarseMax in MAXLEVELS-1? then take all of them, and remove them from the PVS list.
      vNextToSearch = avPVS[mCurrentKF.MAXLEVELS-1];
      avPVS[mCurrentKF.MAXLEVELS-1].clear();
    }
    else
    { // ..otherwise choose nCoarseMax at random, again removing from the PVS list.
      for(unsigned int i=0; i<nCoarseMax; i++)
        vNextToSearch.push_back(avPVS[mCurrentKF.MAXLEVELS-1][i]);
      avPVS[mCurrentKF.MAXLEVELS-1].erase(avPVS[mCurrentKF.MAXLEVELS-1].begin(), avPVS[mCurrentKF.MAXLEVELS-1].begin() + nCoarseMax);
    }

    // If didn't source enough from MAXLEVELS-1, get some from MAXLEVELS-2... same as above.
    if(vNextToSearch.size() < nCoarseMax)
    {
      unsigned int nMoreCoarseNeeded = nCoarseMax - vNextToSearch.size();
      if(avPVS[mCurrentKF.MAXLEVELS-2].size() <= nMoreCoarseNeeded)
      {
        vNextToSearch = avPVS[mCurrentKF.MAXLEVELS-2];
        avPVS[mCurrentKF.MAXLEVELS-2].clear();
      }
      else
      {
        for(unsigned int i=0; i<nMoreCoarseNeeded; i++)
          vNextToSearch.push_back(avPVS[mCurrentKF.MAXLEVELS-2][i]);
        avPVS[mCurrentKF.MAXLEVELS-2].erase(avPVS[mCurrentKF.MAXLEVELS-2].begin(), avPVS[mCurrentKF.MAXLEVELS-2].begin() + nMoreCoarseNeeded);
      }
    }
    // Now go and attempt to find these points in the image!
    unsigned int nFound = SearchForPoints(vNextToSearch, nCoarseRange, *gvnCoarseSubPixIts);

    vIterationSet = vNextToSearch;  // Copy over into the to-be-optimised list.
    if(nFound >= *gvnCoarseMin)  // Were enough found to do any meaningful optimisation?
    {
      mbDidCoarse = true;
      for(int iter = 0; iter<10; iter++) // If so: do ten Gauss-Newton pose updates iterations.
      {
        // incorporate a coarse estimation of the non-rigid model
        // based on the low level points
        if(iter != 0)
        { // Re-project the points on all but the first iteration.
          for(unsigned int i=0; i<vIterationSet.size(); i++)
            if(vIterationSet[i]->bFound)
              vIterationSet[i]->ProjectAndDerivs(mse3CamFromWorld, mCamera);
        }
        for(unsigned int i=0; i<vIterationSet.size(); i++)
          if(vIterationSet[i]->bFound)
            vIterationSet[i]->CalcJacobian(mse3CamFromWorld, mse3CamFromWorldPrev);
        double dOverrideSigma = 0.0;
        // Hack: force the MEstimator to be pretty brutal
        // with outliers beyond the fifth iteration.
        if(iter > 5)
          dOverrideSigma = 1.0;

        // Calculate and apply the pose update...
        Vector<6> v6Update = CalcPoseUpdate(vIterationSet, dOverrideSigma);
        mse3CamFromWorld = SE3<>::exp(v6Update) * mse3CamFromWorld;
      }
    }
  }

  // So, at this stage, we may or may not have done a coarse tracking stage.
  // Now do the fine tracking stage. This needs many more points!

  int nFineRange = nCoarseRange/2; // Pixel search range for the fine stage.
  if(mbDidCoarse)       // Can use a tighter search if the coarse stage was already done.
    nFineRange = 5;

  // What patches shall we use this time? The high-level ones are quite important,
  // so do all of these, with sub-pixel refinement.
  {
    int l = mCurrentKF.MAXLEVELS - 1;
    for(unsigned int i=0; i<avPVS[l].size(); i++)
      avPVS[l][i]->ProjectAndDerivs(mse3CamFromWorld, mCamera);
    SearchForPoints(avPVS[l], nFineRange, 8);
    for(unsigned int i=0; i<avPVS[l].size(); i++)
      // Again, plonk all searched points onto the (maybe already populate) vIterationSet.
      vIterationSet.push_back(avPVS[l][i]);
  }

  // All the others levels: Initially, put all remaining potentially visible patches onto vNextToSearch.
  vNextToSearch.clear();
  for(int l=mCurrentKF.MAXLEVELS - 2; l>=0; l--)
    for(unsigned int i=0; i<avPVS[l].size(); i++)
      vNextToSearch.push_back(avPVS[l][i]);

  // But we haven't got CPU to track _all_ patches in the map - arbitrarily limit
  // ourselves to 1000, and choose these randomly.
  static gvar3<int> gvnMaxPatchesPerFrame("Tracker.MaxPatchesPerFrame", 1000, SILENT);
  int nFinePatchesToUse = *gvnMaxPatchesPerFrame - vIterationSet.size();
  if((int) vNextToSearch.size() > nFinePatchesToUse)
  {
    random_shuffle(vNextToSearch.begin(), vNextToSearch.end());
    vNextToSearch.resize(nFinePatchesToUse); // Chop!
  }

  // If we did a coarse tracking stage: re-project and find derivs of fine points
  if(mbDidCoarse)
    for(unsigned int i=0; i<vNextToSearch.size(); i++)
      vNextToSearch[i]->ProjectAndDerivs(mse3CamFromWorld, mCamera);

  // Find fine points in image:
  SearchForPoints(vNextToSearch, nFineRange, 4);
  // And attach them all to the end of the optimisation-set.
  for(unsigned int i=0; i<vNextToSearch.size(); i++)
    vIterationSet.push_back(vNextToSearch[i]);

  // Again, ten gauss-newton pose update iterations.
  Vector<6> v6LastUpdate;
  v6LastUpdate = Zeros;
  mpInterpolation->Wrap(vIterationSet);
  double dRMSErrorAnt = mdRMSError = computeProjError_RMS(vIterationSet);
  double bestErr=1e6;
  unsigned int nbases = getNBasis();
  SE3<> bestPose;

  double *pdCoeffsCopy = new double[nbases];
  Vector<Dynamic,double,Reference> vBestCoefs(pdCoeffsCopy, nbases);
  double diff = -1;

  string posemethod = GV3::get<string>("posemethod","ptam");
  bool forceRigid = !(GV3::get<string>("forceRigid","").empty());

  for(int iter = 0; iter<10; iter++)
  {
    bool bNonLinearIteration=true; // For a bit of time-saving: don't do full nonlinear
                                   // reprojection at every iteration - it really isn't necessary!
    if(iter == 0 || iter == 4 || iter == 9 || diff < 0)
      bNonLinearIteration = true;   // Even this is probably overkill, the reason we do many
    else                            // iterations is for M-Estimator convergence rather than
      bNonLinearIteration = false;  // linearisation effects.

    // Again, an M-Estimator hack beyond the fifth iteration.
    double dOverrideSigma = 0.0;
    if(iter > 4)
      dOverrideSigma = 160.0;

    if(bNonLinearIteration)
    {
      string coeffsfile = GV3::get<string>("coeffs","");
      setUpdate3DFlags(vIterationSet);
      if((mdRMSError < err_th && coeffsfile.empty() && diff < 0) && (!forceRigid || (forceRigid && mnFrame<1)))
      {
        CalcCoeffUpdate(vIterationSet, dOverrideSigma);
      }

      for(unsigned int i=0; i<vIterationSet.size(); i++)
        if(vIterationSet[i]->bFound)
          vIterationSet[i]->ProjectAndDerivs(mse3CamFromWorld, mCamera);

      for(unsigned int i=0; i<vIterationSet.size(); i++)
        if(vIterationSet[i]->bFound)
          vIterationSet[i]->CalcJacobian(mse3CamFromWorld, mse3CamFromWorldPrev);
    }
    else
    {
      for(unsigned int i=0; i<vIterationSet.size(); i++)
        if(vIterationSet[i]->bFound)
          vIterationSet[i]->LinearUpdate(v6LastUpdate);
    }

    // Calculate and update pose; also store update vector for linear iteration updates.
    Vector<6> v6Update = CalcPoseUpdate(vIterationSet, dOverrideSigma, iter==9);

    //code for debugging purposes
    if (DEBUG){
      cout << "Gauss Newton Iteration " << iter << endl;
      cout << "Update vector" << endl;

      for (int i=0;i<v6Update.size();i++)
        cout << v6Update[i] << " ";
      cout << endl;
      cout << "Camera pose matrix before update" << endl;
      cout << mse3CamFromWorld << endl;
      cout << endl;
    }
    mse3CamFromWorld = SE3<>::exp(v6Update) * mse3CamFromWorld;

    mdRMSError = computeProjError_RMS(vIterationSet);
    if(DEBUG){
      cout << "Camera pose matrix after update" << endl;
      cout << mse3CamFromWorld << " " << endl;
      cout << endl << endl;
    }
    if(bestErr > mdRMSError)
    {
      bestErr = mdRMSError;
      bestPose = mse3CamFromWorld;
      for(unsigned int i = 0; i < nbases; i++)
      {
          vBestCoefs[i] = mpDefCoefs[i];
      }
    }

    diff = dRMSErrorAnt - mdRMSError;
    if(DEBUG)
      cout << endl << "iter " << iter << " err " << mdRMSError << " err ant " << dRMSErrorAnt << " diff " << diff << endl;
    if((diff<5e-5 && diff>0 && bNonLinearIteration) || (diff <5e-6 && diff >0 && !bNonLinearIteration))
      break;
    if(iter !=0 && iter!=4 && iter!=9 && bNonLinearIteration && diff <0)
      break;
    dRMSErrorAnt = mdRMSError;
    v6LastUpdate = v6Update;
  }

  mse3CamFromWorld = bestPose;
  mdRMSError = bestErr;
  for(unsigned int i=0;i<nbases;i++)
  {
    mpDefCoefs[i] = vBestCoefs[i];
  }
  delete [] pdCoeffsCopy;

  mse3CamFromWorldPrev = mse3CamFromWorld;
  //update track window
  for(unsigned int i=0; i<vIterationSet.size(); i++){
    NRTrackerData *TD = vIterationSet[i];
    if(!TD->bFound)
      continue;
    TD->updateDef3D = true;
    TD->Project(mse3CamFromWorld, mCamera);
    TD->v3MapPrev = TD->v3Map;
  }

  if(mbDraw)
  {
    glPointSize(3);
    glEnable(GL_BLEND);
    glEnable(GL_POINT_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_POINTS);
    for(vector<NRTrackerData*>::reverse_iterator it = vIterationSet.rbegin();
        it!= vIterationSet.rend();
        it++)
    {
      if(! (*it)->bFound)
        continue;
      glColor(gavLevelColors[(*it)->nSearchLevel]);
      glVertex((*it)->v2Image);
    }
    glEnd();
    glDisable(GL_BLEND);

  }
  mpInterpolation->compute2DShape();
  mpInterpolation->drawTriangulation();

  // Update the current keyframe with info on what was found in the frame.
  // Strictly speaking this is unnecessary to do every frame, it'll only be
  // needed if the KF gets added to MapMaker. Do it anyway.
  // Export pose to current keyframe:
  mCurrentKF.se3CfromW = mse3CamFromWorld;

  // Record successful measurements. Use the KeyFrame-Measurement struct for this.
  mCurrentKF.mMeasurements.clear();
  for(vector<NRTrackerData*>::iterator it = vIterationSet.begin();
      it!= vIterationSet.end();
      it++)
  {
    if(! (*it)->bFound)
      continue;
    Measurement m;
    m.v2RootPos = (*it)->v2Found;
    m.nLevel = (*it)->nSearchLevel;
    m.bSubPix = (*it)->bDidSubPix;
    mCurrentKF.mMeasurements[& ((*it)->Point)] = m;
  }

  // Finally, find the mean scene depth from tracked features
  {
    double dSum = 0;
    double dSumSq = 0;
    int nNum = 0;
    for(vector<NRTrackerData*>::iterator it = vIterationSet.begin();
        it!= vIterationSet.end();
        it++)
    {
      if((*it)->bFound)
      {
        double z = (*it)->v3Cam[2];
        dSum+= z;
        dSumSq+= z*z;
        nNum++;
      }
    }

    if(nNum > 20)
    {
      mCurrentKF.dSceneDepthMean = dSum/nNum;
      mCurrentKF.dSceneDepthSigma = sqrt((dSumSq / nNum) - (mCurrentKF.dSceneDepthMean) * (mCurrentKF.dSceneDepthMean));
    }
  }

  //clean up old mMapTrails: not existing ones...
  for(unsigned int i=0; i<mMap.vpPointsTrash.size(); i++)
  {
    NRMapPoint *p = mMap.vpPointsTrash[i];
    NRTrackerData* pTD = p->pTData;
    if(mMapTrail.find(pTD)!= mMapTrail.end())
      mMapTrail.erase(pTD);
  }

  mdRMSError = computeBestKFProjError_RMS(vIterationSet);
  mPreviousFrameKF = mCurrentKF;
  if (mbJustRecoveredSoUseCoarse)
    mbJustRecoveredSoUseCoarse = false;
}

/**
 * Instead of the original behaviour, perform trail tracking.
 *
 * @param vTD
 * @param nRange
 * @param nSubPixIts
 * @return
 */
int VisualNRTrackerV1::SearchForPoints(vector<NRTrackerData*> &vTD, int nRange, int nSubPixIts)
{
  int nFound = 0;
  int recovered = 0;

  //as this part of the code will be concurrent mMapTrail must be protected
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for(unsigned int i=0; i<vTD.size(); i++)   // for each point..
  {
    NRTrackerData* TD = vTD[i];

    PatchFinder &Finder = TD->Finder;
    Finder.MakeTemplateCoarseCont(TD->Point);
    TD->bFound = false;
    if(Finder.TemplateBad())
    {
      TD->bInImage = TD->bPotentiallyVisible = TD->bFound = false;
      continue;
    }
    manMeasAttempted[Finder.GetLevel()]++;  // Stats for tracking quality assessment

    bool bFound = Finder.FindPatchCoarse(ir(TD->v2Image), mCurrentKF, nRange);
    TD->bSearched = true;
    if(!bFound)
    {
      //second chance: trail tracking. Even the template could be wrong, but
      //the trails are saved, try to look for them using trail tracking.
      if(!SecondChanceSearch(TD, nRange)){
        TD->bFound = false;
      }else{
        recovered++;
        nFound++;
        manMeasFound[Finder.GetLevel()]++;
        if(TD->nSearchLevel<2)
        {
          Trail &t = mMapTrail[TD];
          ImageRef ir = t.irCurrentPos*LevelScale(TD->nSearchLevel);
          TD->v2Found[0] = ir.x; TD->v2Found[1] = ir.y;
          TD->bFound = true;
        }
      }
      continue;
    }

    TD->bFound = true;
    TD->dSqrtInvNoise = (1.0 / Finder.GetLevelScale());

    nFound++;
    manMeasFound[Finder.GetLevel()]++;

    // Found the patch in coarse search - are Sub-pixel iterations wanted too?
    if(nSubPixIts > 0)
    {
      TD->bDidSubPix = true;
      Finder.MakeSubPixTemplate();
      bool bSubPixConverges = Finder.IterateSubPixToConvergence(mCurrentKF, nSubPixIts);
      if(!bSubPixConverges)
      {
        // If subpix doesn't converge, the patch location is probably very dubious!
        //give it a second chance with miniPatch?
        if(!SecondChanceSearch(TD, nRange))
        {
          TD->bFound = false;
          nFound--;
          manMeasFound[Finder.GetLevel()]--;
        }else{
          recovered++;
          //Transfer the information from the track to the v2Found field
#ifdef _OPENMP
          #pragma omp critical(dataaccess)
#endif
          if(TD->nSearchLevel<2)
          {
            Trail &t = mMapTrail[TD];
            ImageRef ir = t.irCurrentPos*LevelScale(TD->nSearchLevel);
            TD->v2Found[0] = ir.x;TD->v2Found[1] = ir.y;
            TD->bFound = true;
          }
        }
        continue;
      }
      TD->v2Found = Finder.GetSubPixPos();
    }
    else
    {
      TD->v2Found = Finder.GetCoarsePosAsVector();
      TD->bDidSubPix = false;
    }

    //set the trail last position to the one found by the normal tracker
    //if it is not created, this piece of code, creates it
#ifdef _OPENMP
    #pragma omp critical(dataaccess)
#endif
    if(TD->bFound && TD->nSearchLevel<2){ //does it make sense?
      Trail & t = mMapTrail[TD];
      t.irCurrentPos = ir(TD->v2Found)/LevelScale(TD->nSearchLevel);
      //last check-> was the initial position initialized? no-> fill it
      if(t.irInitialPos[0]<=1 && t.irInitialPos[1]<=1)
      {
        t.irInitialPos=t.irCurrentPos;
        t.nSearchLevel=TD->nSearchLevel;
        t.mPatch.SampleFromImage(t.irInitialPos,
                                 mCurrentKF.aLevels[TD->nSearchLevel].im);
      }
      //update search level coordinates of the point, and the patch
      if(t.nSearchLevel!=TD->nSearchLevel)
      {
        t.irInitialPos *= LevelScale(t.nSearchLevel);
        t.irInitialPos /= LevelScale(TD->nSearchLevel);
        t.mPatch.SampleFromImage(t.irCurrentPos,
                                 mCurrentKF.aLevels[TD->nSearchLevel].im);
        t.nSearchLevel = TD->nSearchLevel;
      }
    }
  }
  return nFound;
}

/**
 * This function implements trail tracking based on the MapPoints.
 * @param TD
 * @param nRange
 * @return whether the point is found or not
 */
bool VisualNRTrackerV1::SecondChanceSearch(NRTrackerData* TD, int nRange)
{
  map<NRTrackerData *, Trail>::iterator it;
#ifdef _OPENMP
    #pragma omp critical(dataaccess)
#endif
  {
    it = mMapTrail.find(TD);
  }

  if(it == mMapTrail.end() && TD->nSearchLevel>=2)
    return false;
  nRange /= (TD->nSearchLevel+1);
  nRange = nRange < 1 ? 1: nRange;
  Trail &t = it->second;
  ImageRef irStart = t.irCurrentPos;
  ImageRef irEnd = irStart;

  DeformablePatch finder;
  bool found = finder.FindPatch3(irStart, irEnd, mPreviousFrameKF, mCurrentKF, nRange/2);
  unsigned int lim = nRange/2;
  if(found)
  {
    //married matching
    ImageRef irBackWardsFound = t.irCurrentPos; 
    ImageRef irBackWardsFound2 = irBackWardsFound;

    found = finder.FindPatch3(irBackWardsFound, irBackWardsFound2, mCurrentKF, mPreviousFrameKF, nRange/2);
    unsigned int diff1 = (irBackWardsFound2 - irStart).mag_squared();
    unsigned int diff2 = (irBackWardsFound2 - t.irCurrentPos).mag_squared();

    if(found && diff1 > lim && diff2 > lim)
    {
      found = false;
    }
  }
  if(found)
  {
    t.irCurrentPos = irEnd;
  }
  return found;
}

/**
 * Function that forces the update of the 3D points
 * @param vTD
 * @param value
 */
void VisualNRTrackerV1::setUpdate3DFlags(vector<NRTrackerData*> &vTD, bool value)
{
  //there should be a correspondence between boints and basis
  for(unsigned int i=0;i<vTD.size();i++)
  {
    NRTrackerData *TD = vTD[i];
    if(!TD->bFound)
      continue;
    TD->updateDef3D = true;
  }
}

/**
 * This function calls the helper to estimate the pose based on the set of tracks.
 * @param vTD
 * @param dOverrideSigma
 * @param bMarkOutliers
 */
Vector<6> VisualNRTrackerV1::CalcPoseUpdate(vector<NRTrackerData*> &vTD, double dOverrideSigma, bool bMarkOutliers){
  return mStateEstimator.CalcPoseUpdate(vTD,dOverrideSigma,bMarkOutliers, mnFrame<2);
}

/**
 * This function computes the set of the coefficients to get the deformation 
 * state of the reconstructed object
 * 
 * @param vTD
 * @param dOverrideSigma
 * @param bMarkOutliers
 */
void VisualNRTrackerV1::CalcCoeffUpdate(vector<NRTrackerData*> &vTD,
                                            double dOverrideSigma,
                                            bool bMarkOutliers)
{
  mStateEstimator.CalcCoeffUpdate2(vTD, mpInterpolation,
                                  dOverrideSigma, bMarkOutliers, mnFrame<2, &mMap);
}

/**
 * Restart the shape of the last known correct shape, i.e. the rigid shape
 */
void VisualNRTrackerV1::bootstrapDeformation()
{
  cout << "VisualNRTracker: going back to rigid part for better relocalization" << endl;
  //for all the points in the map, if there is a NRTrackerData available, call bootstrapdef
  Vector<Dynamic,double,Reference> vDefCoefs(mpDefCoefs,mBase.getNBasis());
  vDefCoefs = Zeros;
  for(unsigned int i=0;i<mMap.vpPoints.size();i++)
  {
    NRMapPoint *p = mMap.vpPoints[i];
    NRTrackerData *TD = p->pTData;
    if(!TD)
      continue;
    TD->updateDef3D = true;
  }

  for(unsigned int i=0;i<mMap.vpPointsTrash.size();i++)
  {
    NRMapPoint *p = mMap.vpPointsTrash[i];
    NRTrackerData *TD = p->pTData;
    if(!TD)
      continue;
    TD->updateDef3D = true;
  }
}

/*
 * This function takes the DBases loaded from disk, the rigid shape estimation
 * for the frame and the detected set of points. Then it initializes the Map.
 */
void VisualNRTrackerV1::alignMapFromDisk()
{
  // Are there detected points at this stage? Should be, call this function
  // from the propper place
  // This populates the Candidates list, which is Shi-Tomasi thresholded.
  mCurrentKF.MakeKeyFrame_Rest();
  mCurrentKF.se3CfromW = mse3CamFromWorld;
  mFirstKF = mCurrentKF;
  mFirstKF.se3CfromW = mse3CamFromWorld;
  mFirstKF.MakeKeyFrame_Rest();

  mpInterpolation->compute2DShape();

  //Fill it in now? wait for the tracks to be interpolated and matched?...
  vector<pair<double,ImageRef> > vCornersAndSTScores;
  // Copy candidates into a trivially sortable vector
  // so that we can choose the image corners with max ST score
  for(unsigned int i=0; i<mCurrentKF.aLevels[0].vCandidates.size(); i++)
  {
    Candidate &c = mCurrentKF.aLevels[0].vCandidates[i];
    if(!mCurrentKF.aLevels[0].im.in_image_with_border(c.irLevelPos, MiniPatch::mnHalfPatchSize))
      continue;
    // negative so highest score first in sorted list
    vCornersAndSTScores.push_back(pair<double,ImageRef>(-1.0 * c.dSTScore, c.irLevelPos));
  }
  // Sort according to Shi-Tomasi score
  sort(vCornersAndSTScores.begin(), vCornersAndSTScores.end());
  int nToAdd = GV2.GetInt("MaxInitialTrails", 1000, SILENT);

  vector<NRTrackerData *> v = vector<NRTrackerData *>(nToAdd);
  mMap.vpPoints.clear();

  int roiu2 = GV3::get<unsigned int>("roiu2",0);
  //check the point is inside the roi defined by the four points
  int roiu1 = GV3::get<unsigned int>("roiu1",0);
  int roiv1 = GV3::get<unsigned int>("roiv1",0);
  int roiv2 = GV3::get<unsigned int>("roiv2",0);
  int roiu4 = GV3::get<unsigned int>("roiu4",0);
  int roiv4 = GV3::get<unsigned int>("roiv4",0);

  PatchFinder finder;

  unsigned int idx = 0;

  for(unsigned int i = 0; i<vCornersAndSTScores.size() && nToAdd > 0; i++)
  {
    if(!mCurrentKF.aLevels[0].im.in_image_with_border(vCornersAndSTScores[i].second, MiniPatch::mnHalfPatchSize))
      continue;

    Vector<2> am;
    am[0] = (vCornersAndSTScores[i].second[0]-roiu1);
    am[1] = (vCornersAndSTScores[i].second[1]-roiv1);
    Vector<2> ab;
    ab[0] = roiu2-roiu1;
    ab[1] = roiv2-roiv1;
    Vector<2> ad;
    ad[0] = roiu4-roiu1;
    ad[1] = roiv4-roiv1;

    if(DEBUG)
      cout << " det point " << vCornersAndSTScores[i].second
           << " roi: " << roiu1 << "," << roiv1 << " "
           << roiu2 << "," << roiv2 << " " << roiu4 << "," << roiv4
           << " rect components: " << am << " " << ab << " " << ad << endl;

    if (roiu2 && !(0 < am*ab && am*ab <= ab*ab && 0 < am*ad && am*ad< ad*ad)){
      cout << "rejected as it is not within the limits\n";
      continue;
    }

    //Remember that this is dynamic memory. If something is going wrong, it must be deleted
    NRMapPoint *p = new NRMapPoint();
    p->pPatchSourceKF = &mFirstKF;
    p->nSourceLevel = 0;
    p->v3Normal_NC = makeVector(0,0,-1);
    p->irCenter = vCornersAndSTScores[i].second;
    p->v3Center_NC = unproject(mCamera.UnProject(p->irCenter));
    p->v3OneDownFromCenter_NC = unproject(mCamera.UnProject(p->irCenter + ImageRef(0,1)));
    p->v3OneRightFromCenter_NC = unproject(mCamera.UnProject(p->irCenter + ImageRef(1,0)));
    normalize(p->v3Center_NC);
    normalize(p->v3OneDownFromCenter_NC);
    normalize(p->v3OneRightFromCenter_NC);

    Measurement mFirst;
    mFirst.nLevel = 0;
    mFirst.Source = Measurement::SRC_ROOT;
    //measurement found, substitute with the appropiate data
    mFirst.v2RootPos = vec(vCornersAndSTScores[i].second);
    mFirst.bSubPix = true;
    mFirstKF.mMeasurements[p] = mFirst;

    NRTrackerData *TD = new NRTrackerData(p, &mBase, i, mpDefCoefs);
    //get its 3D pos from the projection and interpolating the mesh on the camera ...
    TD->attached = false;
    TD->bDidSubPix = false;
    TD->bFound = true;
    TD->bInImage = true;
    TD->bOutlier = false;
    TD->bPotentiallyVisible = true;
    TD->bSearched = true;
    //all the tracks are in the first level of the pyramid
    TD->dSqrtInvNoise = 1.0;
    TD->v2Found = vec(vCornersAndSTScores[i].second);
    TD->v2Image = TD->v2Found;

    p->pTData = TD;

    v[idx++] = TD; //already allocated
    mMap.vpPoints.push_back(p);
    nToAdd--;
  }

  mpInterpolation->Wrap(v);

  for( unsigned int i = 0; i < v.size(); i++){
    NRTrackerData* TD = v[i];
    if (!TD) continue;
    cout << " rigid point " << i << " after interpolation " << TD->v3Rigid << endl;
  }

  //set the rigid coefficients to the rigid coefficients.
  Vector<Dynamic,double,Reference> vDefCoefs(mpDefCoefs,getNBasis());
  Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpDefCoefs_rigid,getNBasis());
  vDefCoefs = vDefCoefs_rigid;

  // Now we get the interpolated points, project them. 3D shape is updated.
  // Then, fill in the 3D point information and all the stuff you could get
  for(unsigned int i=0; i<v.size();i++)
  {
    //Take all the current detected points, the current pose and current shape, 
    //then compute their centroids in 2D coords and interpolate...
    NRTrackerData* TD = v[i];
    cout << " rigid point " << i << " after interpolation " << TD->v3Rigid << endl;
    NRMapPoint &p = TD->Point;
    TD->Project(mse3CamFromWorld*mBase.mse3RelativePose, mCamera);

    p.RefreshPixelVectors();

    // Do sub-pixel alignment on the second image
    finder.MakeTemplateCoarseNoWarp(p);
    finder.MakeSubPixTemplate();
    finder.SetSubPixPos(TD->v2Found);
  }
  vDefCoefs = Zeros;
  mMap.vpKeyFrames.push_back(&mFirstKF);
  // set the map to a good state
  mMap.bGood = true;

}
