// -*- c++ -*-
// Author: Sebastian Bronte Palacios. sebastian.bronte@depeca.uah.es
// Date: 3 Feb 2011
// NRTracker.cc
//
// Defines the NRTracker, associated classes and aditional structures needed
//
// This class handles some algorithms to perform tracking based on NRSfM theory.
//

#include "EstimationUtils.h"
#include "OpenGL.h"
#include "MEstimator.h"
#include "ShiTomasi.h"
#include "SmallMatrixOpts.h"
#include "PatchFinder.h"
#include "NRTracker.h"
#include "NRTrackerData.h"

#include <cvd/utility.h>
#include <cvd/gl_helpers.h>
#include <cvd/fast_corner.h>
#include <cvd/vision.h>
#include <TooN/wls.h>
#include <gvars3/instances.h>
#include <gvars3/GStringUtil.h>

#include <fstream>
#include <fcntl.h>

#include<iostream>
#include<string>

#include<sys/time.h>
#include <TooN/SVD.h>

using namespace CVD;
using namespace std;
using namespace GVars3;

//NRTracker constructor
//The same as the parent constructor + read the bases necessary to apply NRSfM.

NRTracker::NRTracker(ImageRef &irVideoSize, const ATANCamera &c, Map &m, DBases &bases)
: Tracker(irVideoSize, c, m, 1), base(bases)
{

    areTracksLoadedFromDisk=false;
    areCoeffsLoadedFromDisk=false;
    mpDefCoefs = NULL;
    mpTracks = NULL;

    string tracksdisk = GV3::get<string>("tracks","");
    if(!tracksdisk.empty())
    {
        areTracksLoadedFromDisk=true;
        cout << "Tracks successfully loaded from file" << endl;
    }

    string coeffsdisk = GV3::get<string>("coeffs","");
    if(!coeffsdisk.empty())
    {
        areCoeffsLoadedFromDisk=true;
        cout << "Coeffs successfully loaded from file" << endl;
    }
    // First a pose estimation is needed, to project everything and know 2D-3D correspondences
    // Quick hack (to keep working on it): I'm god and I know everything xD... at least the
    // initial camera pose.
    // TODO: take this initialization out (it will be set to DBases, fixed)
    // and only compute the relative one on the tracker

    Matrix<3> R;
    R[0] = makeVector(GV3::get<double>("r00",1),GV3::get<double>("r01",0),GV3::get<double>("r02",0));
    R[1] = makeVector(GV3::get<double>("r10",0),GV3::get<double>("r11",1),GV3::get<double>("r12",0));
    R[2] = makeVector(GV3::get<double>("r20",0),GV3::get<double>("r21",0),GV3::get<double>("r22",1));

    Vector<3> T = makeVector(GV3::get<double>("t0",0),GV3::get<double>("t1",0),GV3::get<double>("t2",0));

    mse3CamFromWorld.get_rotation() = R;
    mse3CamFromWorld.get_translation() = T;
    mse3CamFromWorld_default = mse3CamFromWorld;

    // First part of basis initialization
    string basedisk = GV3::get<string>("bases","");

    if(!basedisk.empty())
    {
        normalizedModel = GV3::get<bool>("normalizedModel",false);
        cout << "Normalizing the bases? " << normalizedModel << endl;

        unsigned int nbases = base.getNBasis();
        mpDefCoefs = new double[nbases];
        mpDefCoefs_rigid = new double[nbases];
        //allocating memory for coefficient boostrap
        mpDefCoefs_default = new double[nbases];

        mStateEstimator.setCoefficients(mpDefCoefs, mpDefCoefs_rigid, nbases);
        mStateEstimator.setBases(&bases);
    }

    //second part of basis initialization stuff
    if(!basedisk.empty())
    {
        // Commonsense initialization of the coeficients. Assume the rigid (initially given)
        // shape to be close to the actual one when the program starts.
        unsigned int nbases = base.getNBasis();
        Vector<Dynamic,double,Reference> vDefCoefs(mpDefCoefs,nbases);
        Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpDefCoefs_rigid,nbases);
        vDefCoefs = Zeros;
        vDefCoefs_rigid = Zeros;
 
        // Copy results for coefficient boostrap
        for (unsigned int i=0;i<nbases;i++)
            mpDefCoefs_default[i] = mpDefCoefs[i];
        cout << "Data filled up successfully" << endl;
    }
    applyMotionModel=false;
}

NRTracker::~NRTracker(){
    if(base.getNBasis())
    {
        delete [] mpDefCoefs;
        delete [] mpDefCoefs_default;
        delete [] mpDefCoefs_rigid;
    }
}

void NRTracker::Reset()
{
    Tracker::Reset();
    //leave the coefficients as they were
    for(unsigned int i=0; i<base.getNBasis();i++)
    {
      mpDefCoefs[i] = mpDefCoefs_default[i];
    }
}

/**
 * NRTracker TrackFrame functions
 * This overloaded function selects the destination function, whether call the
 * parent function or the one in this class
 */
void NRTracker::TrackFrame(Image<byte> &imFrame, bool bDraw)
{
    mbDraw = bDraw;
    mMessageForUser.str("");   // Wipe the user message clean

    static gvar3<int> gvnUseSBI("Tracker.UseRotationEstimator", 0, SILENT); //not using sbi rotation stimation until everything is working as it was doing previously
    mbUseSBIInit = *gvnUseSBI;

    if (areTracksLoadedFromDisk)
    {
      TrackFrameFromGivenTracks();
    }
    else
      Tracker::TrackFrame(imFrame,bDraw);

    while(!mvQueuedCommands.empty())
    {
      GUICommandHandler(mvQueuedCommands.begin()->sCommand, mvQueuedCommands.begin()->sParams);
      mvQueuedCommands.erase(mvQueuedCommands.begin());
    }
}

/**
 * TrackFrame function version for text-based input data.
 * @param tracks: the aligned tracking information of the 2D tracking points over time.
 * @param coeffs: the coefficient vector for the whole sequence
 *                (if loaded with the corresponding option).
 * @param bDraw: whether to draw the results for the current frame or not.
 */
void NRTracker::TrackFrame(vector<Vector<2> > &tracks, vector<bool> &visibility, vector<double> &coeffs, bool bDraw)
{
    mbDraw = bDraw;
    mpTracks = &tracks;

    if(mnFrame==0)
      alignMapFromDisk();

    //asign visibility for tracks
    if(!visibility.empty())
    {
        manMeasFound[0]=0;
        //manMeasAttempted[0] = mMap.vpPoints.size();
        for(unsigned int i=0; i<mMap.vpPoints.size(); i++)
        {
            NRMapPoint *p = mMap.vpPoints[i];
            bool tmpvis = visibility[i];
            //For text based data assume just 1 level of the pyramid
            //if(tmpvis) manMeasFound[0]++;
            p->pTData->bFound = tmpvis;
        }
    }

    if(!coeffs.empty())
    {
        for(unsigned int i=0;i<base.getNBasis();i++)
        {
            mpDefCoefs[i] = coeffs[i];
        }
    }

    mMessageForUser.str("");   // Wipe the user message clean

    static gvar3<int> gvnUseSBI("Tracker.UseRotationEstimator", 0, SILENT); //not using sbi rotation stimation until everything is working as it was doing previously
    mbUseSBIInit = *gvnUseSBI;

    TrackFrameFromGivenTracks();

    while(!mvQueuedCommands.empty())
    {
      GUICommandHandler(mvQueuedCommands.begin()->sCommand, mvQueuedCommands.begin()->sParams);
      mvQueuedCommands.erase(mvQueuedCommands.begin());
    }
}

void NRTracker::TrackFrameFromGivenTracks()
{
    if(mMap.IsGood())
    {
        if(mnLostFrames < 3)  // .. but only if we're not lost!
        {
            if(normalizedModel)
            {
                UpdateDeformationModelSize();
            }
            if(applyMotionModel)
            {
                if(mbUseSBIInit)
                    CalcSBIRotation();
                ApplyMotionModel();
            }

            TrackMap();               //  This line does the main tracking work.

            if(applyMotionModel)
                UpdateMotionModel();
            // Check if we're lost or if tracking is poor.
            // It depends on the tracking type we are considering.
            Tracker::AssessTrackingQuality();  
        }
        else
        {
            mMap.bGood=false;
        }
    }
    else
    {
        mMessageForUser << "** Attempting recovery **.";
        if(AttemptRecovery())
        {
          TrackMap();
          Tracker::AssessTrackingQuality();
        }
        if(mTrackingQuality==GOOD)
        {
            mMap.bGood=true;
            mnLostFrames=0;
        }
        else mMap.bGood=false;
    }

    { // Provide some feedback for the user:
        mMessageForUser << "Tracking Map, quality ";
        if(mTrackingQuality == GOOD)  mMessageForUser << "good.";
        if(mTrackingQuality == DODGY) mMessageForUser << "poor.";
        if(mTrackingQuality == BAD)   mMessageForUser << "bad.";
        mMessageForUser << " Found:";
        for(int i=0; i<mCurrentKF.MAXLEVELS; i++) mMessageForUser << " " << manMeasFound[i] << "/" << manMeasAttempted[i];
        mMessageForUser << " Map: " << mMap.vpPoints.size() << "P, " << mMap.vpKeyFrames.size() << "KF, " << base.getNBasis() << "B";
    }

    mnFrame++;
}

/**
 * called if tracking is lost. As easy as give the last good pose to the system
 * or the first one if it's near the begining
 */
bool NRTracker::AttemptRecovery()
{
    //copy default pose
    mse3CamFromWorld = mse3CamFromWorld_default;

    //copy default coeficients
    for(unsigned int i=0;i<base.getNBasis();i++)
        mpDefCoefs[i]=mpDefCoefs_default[i];

    return true;
}

/**
 * Funtion to align the map if it is loaded from disk.
 * this function performs data association, fill as many members as possible
 * to avoid crashes
 */
void NRTracker::alignMapFromDisk()
{
    //the corners have already been filled in, so the only thing to do keep
    //going and using the info contained in the current keyframe
    KeyFrame *pkFirst = new KeyFrame();
    pkFirst->bFixed = true;
    pkFirst->se3CfromW = SE3<>();

    //get the tracks for the first frame
    vector<Vector<2> > vTracks = *mpTracks;

    for(unsigned int i=0;i<base.getNPoints();i++)
    {
        //fill in all the points using whether the 1st basis or the average one.
        //data association is direct since the same index is applied to 2D and 3D data
        //checking if there is some da1ta already loaded (in the case of loading
        //both rigid average shape and basis, the rigid one is loaded first)
        NRMapPoint *p;
        p = new NRMapPoint();
        p->v3WorldPos = Zeros;
        for(unsigned int k=0;k<base.getNBasis();k++)
          p->v3WorldPos += base.getBasisPoint(k,i)*mpDefCoefs[k];
        p->nSourceLevel = 0;
        mMap.vpPoints.push_back(p);
    }
    //mMapMaker.fillIn3DData(base, mpDefCoefs);

    Matrix<> & rigidSet = *(base.mmRigidShape);

    for(unsigned int i=0; i<mMap.vpPoints.size(); i++)
    {
        NRMapPoint *p = mMap.vpPoints[i];

        //3D point already filled up, then the rest of the members
        p->pPatchSourceKF = pkFirst;
        p->nSourceLevel = 0;
        p->v3Normal_NC = makeVector( 0,0,-1);

        // Ensure that this map point has an associated NRTrackerData struct.
        if(!p->pTData) p->pTData = new NRTrackerData(p,&base,i,mpDefCoefs);
        NRTrackerData *TData = p->pTData;
        if (base.isRigidFromFile()) TData->v3Rigid = rigidSet[i];
        //copying for the next time the value of the non rigid term.
        TData->v3Def = p->v3WorldPos;

        TData->updateDef3D = false;
        //compute the whole 3D position for the next stages
        p->v3WorldPos += TData->v3Rigid;

        // Project according to current view
        TData->Project(mse3CamFromWorld, mCamera);
        if(!TData->bInImage)
            continue;

        //point projected derivs obtained, find it on the image, should be the same. 
        //anyway is not needed, no patches are necessary, no warping matrix will be needed
        TData->nSearchLevel = 0;
        TData->dSqrtInvNoise = 1.0;

        TData->bSearched = true;
        TData->bFound = true;

        //fill in the rest of the necesary data
        p->irCenter = ir(vTracks[i]);
        p->v3Center_NC = unproject(mCamera.UnProject(p->irCenter));
        p->v3OneRightFromCenter_NC = unproject(mCamera.UnProject(p->irCenter + ImageRef(1,0)));
        normalize(p->v3Center_NC);
        normalize(p->v3OneDownFromCenter_NC);
        normalize(p->v3OneRightFromCenter_NC);

        Measurement mFirst;
        mFirst.nLevel = 0;
        mFirst.Source = Measurement::SRC_ROOT;
	//measurement found, substitute with the appropiate data
        mFirst.v2RootPos = vTracks[i];
        mFirst.bSubPix = true;
        pkFirst->mMeasurements[p] = mFirst;
    }

    //adding the keyframe to the list of keyframes
    mMap.vpKeyFrames.push_back(pkFirst);

    //let the Tracker thread know that the map is good
    mMap.bGood = true;
}

/**
 * When a map is available, performs the necessary operations to track the 3D points
 * with the detected ones in 2D
 */
void NRTracker::TrackMap()
{
    // measuring computation time
    struct timeval t1,t2;

    // Some accounting which will be used for tracking quality assessment:
    for(int i=0; i<mCurrentKF.MAXLEVELS; i++)
        manMeasAttempted[i] = manMeasFound[i] = 0;

    // All the data is given but the tracks are loaded sequencially
    vector<Vector<2> > tracksInFrame = *mpTracks;

    // The Potentially-Visible-Set (PVS) is split into pyramid levels.
    vector<NRTrackerData*> avPVS[mCurrentKF.MAXLEVELS];
    for(int i=0; i<mCurrentKF.MAXLEVELS; i++)
        avPVS[i].reserve(500);

    gettimeofday(&t1,NULL);

    // For all points in the map..
    for(unsigned int i=0; i<mMap.vpPoints.size(); i++)
    {
        NRMapPoint &p= *(mMap.vpPoints[i]);
        // Ensure that this map point has an associated NRTrackerData struct.
        if(!p.pTData) p.pTData = new NRTrackerData(&p, &base, i,mpDefCoefs);
        NRTrackerData *TData = p.pTData;

        // Project according to current view, and if it is not within the image, skip.
        TData->Project(mse3CamFromWorld, mCamera);
        if(!TData->bInImage)
            continue;
        TData->v2Found = tracksInFrame[i];
        // Calculate camera projection derivatives of this point. Should be convenient to compute the things related to NR stuff here
        TData->GetDerivsUnsafe(mCamera);

        // And check what the PatchFinder (included in NRTrackerData) makes of the mappoint in this view.. -> there is no appearance related content, so, let's put this appart.
        TData->nSearchLevel = 0; //TData.Finder.CalcSearchLevelAndWarpMatrix(TData.Point, mse3CamFromWorld, TData.m2CamDerivs); //when feature detection is enabled, just enable this line again.
        TData->dSqrtInvNoise = 1.0;

        TData->v2Error_CovScaled = TData->dSqrtInvNoise*(TData->v2Found-TData->v2Image);
        // Otherwise, this point is suitable to be searched in the current image! Add to the PVS.
        TData->bSearched = true;
        //TData->bFound = true;
        avPVS[TData->nSearchLevel].push_back(TData);
        manMeasAttempted[0]++; //adding one to the measAttempted array
        if(TData->bFound)
            manMeasFound[0]++;
    };

    gettimeofday(&t2,NULL);
    if(DEBUG)
        cout << "map points projection time: " << (t2.tv_sec-t1.tv_sec)*1000 + (t2.tv_usec-t1.tv_usec)*0.001 << " ms" << endl;

    // The next two data structs contain the list of points which will be next
    // searched for in the image, and then used in pose update.
    vector<NRTrackerData*> vNextToSearch;
    vIterationSetNR.clear();

    // Tunable parameters to do with the coarse tracking stage:
    static gvar3<unsigned int> gvnCoarseMin("NRTracker.CoarseMin", 20, SILENT);   // Min number of large-scale features for coarse stage
    static gvar3<unsigned int> gvnCoarseMax("NRTracker.CoarseMax", 60, SILENT);   // Max number of large-scale features for coarse stage
    static gvar3<unsigned int> gvnCoarseRange("NRTracker.CoarseRange", 30, SILENT);       // Pixel search radius for coarse features
    static gvar3<int> gvnCoarseSubPixIts("NRTracker.CoarseSubPixIts", 8, SILENT); // Max sub-pixel iterations for coarse features
    static gvar3<int> gvnCoarseDisabled("NRTracker.DisableCoarse", 0, SILENT);    // Set this to 1 to disable coarse stage (except after recovery)
    static gvar3<double> gvdCoarseMinVel("NRTracker.CoarseMinVelocity", 0.006, SILENT);  // Speed above which coarse stage is used.

    unsigned int nCoarseMax = *gvnCoarseMax;
    unsigned int nCoarseRange = *gvnCoarseRange;

    mbDidCoarse = false;

    // Set of heuristics to check if we should do a coarse tracking stage.
    if(*gvnCoarseDisabled ||
      mdMSDScaledVelocityMagnitude < *gvdCoarseMinVel  ||
      nCoarseMax == 0)
    if(mbJustRecoveredSoUseCoarse)
    {
        nCoarseMax *=2;
        nCoarseRange *=2;
        mbJustRecoveredSoUseCoarse = false;
    };

    {
        int l = mCurrentKF.MAXLEVELS - 1;
        for(unsigned int i=0; i<avPVS[l].size(); i++)
            vIterationSetNR.push_back(avPVS[l][i]);  // Again, plonk all searched points onto the (maybe already populate) vIterationSetNR.
    };

    // All the others levels: Initially, put all remaining potentially visible patches onto vNextToSearch.
    vNextToSearch.clear();

    // But we haven't got CPU to track _all_ patches in the map - arbitrarily limit
    // ourselves to 1000, and choose these randomly.
    static gvar3<int> gvnMaxPatchesPerFrame("NRTracker.MaxPatchesPerFrame", 1000, SILENT);
    int nFinePatchesToUse = *gvnMaxPatchesPerFrame - vIterationSetNR.size();
    if((int) vNextToSearch.size() > nFinePatchesToUse)
    {
        random_shuffle(vNextToSearch.begin(), vNextToSearch.end());
        vNextToSearch.resize(nFinePatchesToUse); // Chop!
    };

    // If we did a coarse tracking stage: re-project and find derivs of fine points
    if(mbDidCoarse)
        for(unsigned int i=0; i<vNextToSearch.size(); i++)
            vNextToSearch[i]->ProjectAndDerivs(mse3CamFromWorld, mCamera);

    Vector<6> vLastUpdate = Zeros;
    Vector<Dynamic,double,Reference> vDefCoefs(mpDefCoefs,base.getNBasis());

    static gvar3<double> gvdErrTh("NRTracker.errorThreshold", 0.3, SILENT);
    double err_th = *gvdErrTh;
    double mean_error_final = 1e5;

    if(DEBUG)
        cout << "Frame number: " << mnFrame << endl;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    string coeffsfile = GV3::get<string>("coeffs","");
    bool forceRigid = !(GV3::get<string>("forceRigid","").empty());

    if(!coeffsfile.empty())
      setUpdate3DFlags(vIterationSetNR);

    // Ten gauss-newton pose update iterations.
    for(int iter = 0; iter<10; iter++)
    {
        bool bNonLinearIteration; // For a bit of time-saving: don't do full nonlinear
                                  // reprojection at every iteration - it really isn't necessary!
        double mean_error = 1e5;

        //first let's get the non-linear update working properly, and then... we'll go for the linear update
        if(iter == 0 || iter == 4 || iter == 9)
            bNonLinearIteration = true;   // Even this is probably overkill, the reason we do many
        else                            // iterations is for M-Estimator convergence rather than
            bNonLinearIteration = false;  // linearisation effects.

        if(bNonLinearIteration)
        {
            mean_error = computeProjError_RMS(vIterationSetNR);
            if(DEBUG)
                cout << "mean projection error before step 1: " << mean_error << (mean_error>err_th ? " not converged" : " converged") << endl;

            //As most of the papers do, update the coefficients independently
            //the rotation and translation estimation is good enought to start the update
            gettimeofday(&t1,NULL);

            if(mean_error<err_th && coeffsfile.empty() && (!forceRigid || (forceRigid && mnFrame<1)))
            {
                CalcCoeffUpdate(vIterationSetNR);
            }

            gettimeofday(&t2,NULL);
            if(DEBUG)
                cout << "Computing coeffs time: " << (t2.tv_sec-t1.tv_sec)*1000 + (t2.tv_usec-t1.tv_usec)*0.001 << " ms" << endl;

            gettimeofday(&t1,NULL);
            //update first the coordinates for each point in an exhaustive way
            for(unsigned int i=0; i<vIterationSetNR.size(); i++)
                if(vIterationSetNR[i]->bFound)
                    vIterationSetNR[i]->ProjectAndDerivs(mse3CamFromWorld, mCamera);
        }
        else
        {
            gettimeofday(&t1,NULL);
            //Update first the coordinates for each point
            for(unsigned int i=0; i<vIterationSetNR.size(); i++)
                if(vIterationSetNR[i]->bFound)
                    vIterationSetNR[i]->LinearUpdate(vLastUpdate);
        };

        gettimeofday(&t2,NULL);
        if(DEBUG)
            cout << "Point projection time: " << (t2.tv_sec-t1.tv_sec)*1000 + (t2.tv_usec-t1.tv_usec)*0.001 << " ms" << endl;

        gettimeofday(&t1,NULL);

        //try to change the threshold by something estimated dynamically
        if(DEBUG)
            cout << "mean projection error before step 2: " << mean_error << (mean_error>err_th ? " not converged" : " converged, refining") << endl;
        mean_error_final = mean_error;

        //here you compute the Jacobian por each point
        if(bNonLinearIteration)
            for(unsigned int i=0; i<vIterationSetNR.size(); i++)
                if(vIterationSetNR[i]->bFound)
                    vIterationSetNR[i]->CalcJacobian();

        // Again, an M-Estimator hack beyond the fifth iteration.
        double dOverrideSigma = 0.0;
        if(iter > 5)
            dOverrideSigma = 16.0;

        // Calculate and update pose; also store update vector for linear iteration updates.
        Vector<> vUpdate =
                CalcPoseUpdate(vIterationSetNR, dOverrideSigma, false/*iter==9*/);
        if(DEBUG)
            cout << "Gauss Newton Iteration " << iter << endl;

        //code for debugging purposes
        if(DEBUG)
        {
            cout << "Update vector" << vUpdate << endl;

            cout << "Camera pose matrix before update" << endl;
            cout << mse3CamFromWorld << endl << endl;

            cout << "Current frame deformation coefficients" << endl;
        }
        Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpDefCoefs_rigid,base.getNBasis());
        if(DEBUG)
            cout << vDefCoefs+vDefCoefs_rigid << endl << endl;

        //update the pose state
        mse3CamFromWorld = SE3<>::exp(vUpdate) * mse3CamFromWorld;

        if(DEBUG)
        {
            cout << "Camera pose matrix after update" << endl;
            cout << mse3CamFromWorld << endl << endl;
        }
        vLastUpdate = vUpdate;

        gettimeofday(&t2,NULL);
        if(DEBUG)
            cout << "rest of gauss newton iteration time: " << (t2.tv_sec-t1.tv_sec)*1000 + (t2.tv_usec-t1.tv_usec)*0.001 << " ms" << endl;

        //code to visualize the update process performed by this Gauss-newton steps.
        if(mbDraw)
        {
            //configure output
            glDrawPixels(mCurrentKF.aLevels[0].im);
            glPointSize(3);
            glEnable(GL_BLEND);
            glEnable(GL_POINT_SMOOTH);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glBegin(GL_POINTS);
            glColor(makeVector(0.0, 0.2+0.08*iter, 0.5+0.05*iter));

            for(unsigned int i=0;i<vIterationSetNR.size();i++)
            {
                if(vIterationSetNR[i]->bFound)
                    glVertex(vIterationSetNR[i]->v2Image);
            }

            glEnd();
            glDisable(GL_BLEND);

        }
    }

    if(mean_error_final < err_th)
        applyMotionModel=false;
    else
        applyMotionModel=true;
    cout << "applyMotionModel " << applyMotionModel << endl;

    //showing results
    if(mbDraw)
    {
        glDrawPixels(mCurrentKF.aLevels[0].im);
        glPointSize(6);
        glEnable(GL_BLEND);
        glEnable(GL_POINT_SMOOTH);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBegin(GL_POINTS);
        for(vector<NRTrackerData*>::reverse_iterator it = vIterationSetNR.rbegin();
            it!= vIterationSetNR.rend();
            it++)
        {
            //print predictions / model visualization
            glColor(gavLevelColors[(*it)->nSearchLevel]);
            glVertex((*it)->v2Image);

            if(! (*it)->bFound)
                continue;

            //print given tracks
            glColor(makeVector(1.0, 1.0, 1.0));
            glVertex((*it)->v2Found);
        }

        glEnd();
        glDisable(GL_BLEND);

    }

    // Update the current keyframe with info on what was found in the frame.
    // Strictly speaking this is unnecessary to do every frame, it'll only be
    // needed if the KF gets added to MapMaker. Do it anyway.
    // Export pose to current keyframe:
    mCurrentKF.se3CfromW = mse3CamFromWorld;

    // Record successful measurements. Use the KeyFrame-Measurement struct for this.
    mCurrentKF.mMeasurements.clear();
    for(vector<NRTrackerData*>::iterator it = vIterationSetNR.begin();
        it!= vIterationSetNR.end();
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

        for(vector<NRTrackerData*>::iterator it = vIterationSetNR.begin();
            it!= vIterationSetNR.end();
            it++)
        {
            if((*it)->bFound)
            {
                double z = (*it)->v3Cam[2];
                dSum+= z;
                dSumSq+= z*z;
                nNum++;
            };
        }

        if(nNum > 20)
        {
            mCurrentKF.dSceneDepthMean = dSum/nNum;
            mCurrentKF.dSceneDepthSigma = sqrt((dSumSq / nNum) - (mCurrentKF.dSceneDepthMean) * (mCurrentKF.dSceneDepthMean));
        }
    }
}

//CalCoeffUpdate function
//It performs coefficient update based on the given tracks, the already known
//basis and the previous pose estimated.
void NRTracker::CalcCoeffUpdate(vector<NRTrackerData*> &vTD)
{
    mStateEstimator.CalcCoeffUpdate(vTD);
    setUpdate3DFlags(vTD);
}

void NRTracker::setUpdate3DFlags(vector<NRTrackerData*> &vTD)
{
    for(unsigned int i=0;i<vTD.size();i++) //there should be a correspondence between boints and basis
    {
        NRTrackerData *TD = vTD[i];
        TD->updateDef3D = true;
    }
}

void NRTracker::rigidShapeCoeffs()
{
    //solve the problem using least squares given all the 3D correspondences
    //Create matrix A from the already loaded bases.
    mStateEstimator.rigidShapeCoeffs();
    normalizeRigidCoefficients();
}

void NRTracker::normalizeCoefficients()
{
    if(normalizedModel)
    {
        double sum=0;
        for(unsigned int i=0;i<base.getNBasis();i++)
        {
            mpDefCoefs[i] = mpDefCoefs[i]*base.getComponentEnergy(i);
            sum += mpDefCoefs[i];
        }
        currentShapeEnergy = sum;
        ApplyDeformationModelSize();
        for(unsigned int i=0;i<base.getNBasis();i++)
        {
            mpDefCoefs[i] /=currentShapeEnergy;
        }
    }
}

void NRTracker::normalizeRigidCoefficients()
{
    if(normalizedModel)
    {
        //Compute normalization constant
        double sum = 0;
        for(unsigned int i=0;i<base.getNBasis();i++)
        {
            mpDefCoefs_rigid[i] = mpDefCoefs_rigid[i]*base.getComponentEnergy(i);
            sum += mpDefCoefs_rigid[i];
        }
        rigidShapeEnergy = sum;
        currentShapeEnergy = rigidShapeEnergy;
        sum = 1/sum;
        for(unsigned int i=0;i<base.getNBasis();i++)
        {
            mpDefCoefs_rigid[i] *=sum;
        }
    }
}

void NRTracker::ApplyDeformationModelSize()
{
    currentShapeEnergy = rigidShapeEnergy*defModelUpdate;
}

void NRTracker::UpdateDeformationModelSize()
{
    //Assume the size of the object is fixed
    static gvar3<double> gvdbeta("NRTracker.DeformationModelWeight", 0, SILENT);
    double beta = *gvdbeta;
    //this data needs to be normalized from 0 to 1
    defModelUpdate = 1+beta*(currentShapeEnergy/rigidShapeEnergy-1);
}
