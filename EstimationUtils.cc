/**
 * EstimationUtils.cc
 * 
 * Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
 * 
 * Date: 13/01/2013
 * 
 * This file implements EstimationUtils class. It contains functions to perform
 * state computation: pose and coefficient estimation.
 * 
 */

#include "EstimationUtils.h"
#include "NRTrackerData.h"
#include "MEstimator.h"
#include "DBases.h"
#include "Interpolation.h"

#include <TooN/wls.h>
#include <TooN/SVD.h>

#include <gvars3/instances.h>

#include "opencv2/core/core.hpp"
#include "opencv2/calib3d/calib3d.hpp"

#include <vector>

using namespace std;
using namespace TooN;
using namespace GVars3;

/**
 * Class constructor. Iniciates the reference to the needed data
 *
 * @param camera Camera needed for projections
 * @param pose Pose of the current camera
 * @param pCoeffs The pointer to the set of current coefficients
 * @param pRigidCoeffs The pointer to the set of current rigid coefficients
 */
EstimationUtils::EstimationUtils(ATANCamera & camera, SE3<> & pose):
                                 mCamera(camera),mse3Pose(pose),
                                 mpdCoeffs(NULL),mpdRigidCoeffs(NULL),
                                 mnBasis(0), mpBases(NULL)
{
}

/**
 * This function sets the pointers of the coefficients, for the rigid 
 * and non-rigid sets
 * 
 * @param pCoeffs The pointer to the set of current coefficients
 * @param pRigidCoeffs The pointer to the set of current rigid coefficients.
 * @param nBasis The number of basis shapes used.
 */
void EstimationUtils::setCoefficients(double * pCoeffs,
                                      double * pRigidCoeffs,
                                      unsigned int nBasis)
{
  mpdCoeffs = pCoeffs;
  mpdRigidCoeffs = pRigidCoeffs;
  mnBasis = nBasis;
}

/**
 * This function sets the bases for non-rigid computations if needed
 * 
 * @param bases The pointer to the DBases object containing the basis shapes.
 */
void EstimationUtils::setBases(DBases * bases)
{
  mpBases = bases;
}

void EstimationUtils::setPose(SE3<> & se3CamFromPose)
{
  mse3Pose = se3CamFromPose;
}

/**
 * This function estimates the pose update from the mapPoint tracking 
 * information and the sigma limit to include them or not in the results
 * 
 * @param vTD The vector of Tracking Data pointers to get the Tracking data from.
 * @param dOverrideSigma The sigma from which the one from the M-estimator is overriden.
 * @param bMarkOutliers Whether the outlier points should be marked or not.
 * @return pose update vector
 */
Vector<6> EstimationUtils::CalcPoseUpdate(std::vector<NRTrackerData *> &vTD,
                                          double dOverrideSigma,
                                          bool bMarkOutliers,
                                          bool bFirstFrame)
{
  // Which M-estimator are we using?
  int nEstimator = 0;
  static gvar3<string> gvsEstimator("TrackerMEstimator", "Tukey", SILENT);
  if(*gvsEstimator == "Tukey")
    nEstimator = 0;
  else if(*gvsEstimator == "Cauchy")
    nEstimator = 1;
  else if(*gvsEstimator == "Huber")
    nEstimator = 2;
  else
  {
    cout << "Invalid TrackerMEstimator, choices are Tukey, Cauchy, Huber" << endl;
    nEstimator = 0;
    *gvsEstimator = "Tukey";
  }

  // Uncomment to apply the same time smoothing prior factor to the pose.
  //bool tsmoothprior = GV3::get<int>("timesmoothness","0")==1 && !bFirstFrame;
  //double tsmoothFactor = GV3::get<double>("tsmoothfactor","1.0");
  
  // Find the covariance-scaled reprojection error for each measurement.
  // Also, store the square of these quantities for M-Estimator sigma squared estimation.
  vector<double> vdErrorSquared;
  for(unsigned int f=0; f<vTD.size(); f++)
  {
    NRTrackerData &TD = *vTD[f];
    if(!TD.bFound)
      continue;
    //should the v2Error_CovScaled be computed here or previously computed?
    TD.v2Error_CovScaled = TD.dSqrtInvNoise* (TD.v2Found - TD.v2Image);
    vdErrorSquared.push_back(TD.v2Error_CovScaled * TD.v2Error_CovScaled);
  }

  // No valid measurements? Return null update.
  if(vdErrorSquared.size() == 0)
    return makeVector( 0,0,0,0,0,0);

  // What is the distribution of errors?
  double dSigmaSquared;
  if(dOverrideSigma > 0)
    dSigmaSquared = dOverrideSigma; // Bit of a waste having stored the vector of square errors in this case!
  else
  {
    if (nEstimator == 0)
      dSigmaSquared = Tukey::FindSigmaSquared(vdErrorSquared);
    else if(nEstimator == 1)
      dSigmaSquared = Cauchy::FindSigmaSquared(vdErrorSquared);
    else
      dSigmaSquared = Huber::FindSigmaSquared(vdErrorSquared);
  }

  // The TooN WLSCholesky class handles reweighted least squares.
  // It just needs errors and jacobians.
  WLS<6,double,SQSVD> wls;
  wls.add_prior(100.0); // Stabilising prior
  for(unsigned int f=0; f<vTD.size(); f++)
  {
    NRTrackerData &TD = *vTD[f];
    if(!TD.bFound)
      continue;
    Vector<2> &v2 = TD.v2Error_CovScaled;
    double dErrorSq = v2 * v2;
    double dWeight;

    if(nEstimator == 0)
      dWeight= Tukey::Weight(dErrorSq, dSigmaSquared);
    else if(nEstimator == 1)
      dWeight= Cauchy::Weight(dErrorSq, dSigmaSquared);
    else
      dWeight= Huber::Weight(dErrorSq, dSigmaSquared);

    // Inlier/outlier accounting, only really works for cut-off estimators such as Tukey.
    if(dWeight < 0.5)
    {
      if(bMarkOutliers)
        TD.Point.nMEstimatorOutlierCount++;
      TD.setOutlier(true);
      continue;
    }
    else
    {
      if(bMarkOutliers)
        TD.Point.nMEstimatorInlierCount++;
        TD.setOutlier(false);
    }
    Matrix<2,6> &m26Jac = TD.m26Jacobian;
    wls.add_mJ(v2[0], TD.dSqrtInvNoise * m26Jac[0], dWeight); // These two lines are currently
    wls.add_mJ(v2[1], TD.dSqrtInvNoise * m26Jac[1], dWeight); // the slowest bit of poseits

    // Uncomment to apply the same time smoothing prior to the pose.
    //if (tsmoothprior){
    //  wls.add_mJ(0.0, TD.dSqrtInvNoise*TD.priorJacobian, tsmoothFactor*dWeight); //prior on the camera motion
    //}
  }

  wls.compute();
  return wls.get_mu();
}

/**
 * This function is not used, although it is left in case tests with this pose
 * estimation algorithm want to be performed.
 * @param vTD. vector of NRTrackerData pointers.
 * @param dOverrideSigma sigma value from which the estimator is overriden.
 * @param bMarkOutliers whether the outliers should be marked or not.
 * @param bUsePrevPose In case the previous poses wanted to be used to support
 *                     and smooth the current estimation.
 */
SE3<> EstimationUtils::CalcPoseEPnP(std::vector<NRTrackerData *> &vTD,
                                      double dOverrideSigma,
                                      bool bMarkOutliers,
                                      bool bUsePrevPose)
{
  // Maintain the m-estimator to reject points which weights are lower than a threshold
  // adapt the selected points to vector<Point3f>
  // compute the rotation and translation vectors
  // convert from opencv to toon and save the result.
  // modify the corresponding stuff on VisualNRTracker
  // Which M-estimator are we using?
  int nEstimator = 0;
  static gvar3<string> gvsEstimator("TrackerMEstimator", "Tukey", SILENT);
  if(*gvsEstimator == "Tukey")
    nEstimator = 0;
  else if(*gvsEstimator == "Cauchy")
    nEstimator = 1;
  else if(*gvsEstimator == "Huber")
    nEstimator = 2;
  else
  {
    cout << "Invalid TrackerMEstimator, choices are Tukey, Cauchy, Huber" << endl;
    nEstimator = 0;
    *gvsEstimator = "Tukey";
  }

  // Find the covariance-scaled reprojection error for each measurement.
  // Also, store the square of these quantities for M-Estimator sigma squared estimation.
  vector<double> vdErrorSquared;
  for(unsigned int f=0; f<vTD.size(); f++)
  {
    NRTrackerData &TD = *vTD[f];
    if(!TD.bFound)
      continue;
    //should the v2Error_CovScaled be computed here or previously computed?
    TD.v2Error_CovScaled = TD.dSqrtInvNoise* (TD.v2Found - TD.v2Image);
    vdErrorSquared.push_back(TD.v2Error_CovScaled * TD.v2Error_CovScaled);
  }

  // No valid measurements? Return null update.
  if(vdErrorSquared.size() == 0)
    return makeVector(0,0,0,0,0,0);

  // What is the distribution of errors?
  double dSigmaSquared;
  if(dOverrideSigma > 0)
    // Bit of a waste having stored the vector of square errors in this case!
    dSigmaSquared = dOverrideSigma;
  else
  {
    if (nEstimator == 0)
      dSigmaSquared = Tukey::FindSigmaSquared(vdErrorSquared);
    else if(nEstimator == 1)
      dSigmaSquared = Cauchy::FindSigmaSquared(vdErrorSquared);
    else
      dSigmaSquared = Huber::FindSigmaSquared(vdErrorSquared);
  }

  vector<cv::Point2f> projPoints;
  vector<cv::Point3f> modelPoints;
  for(unsigned int f=0; f<vTD.size(); f++)
  {
    NRTrackerData &TD = *vTD[f];
    if(!TD.bFound)
      continue;
    Vector<2> &v2 = TD.v2Error_CovScaled;
    double dErrorSq = v2 * v2;
    double dWeight;

    if(nEstimator == 0)
      dWeight= Tukey::Weight(dErrorSq, dSigmaSquared);
    else if(nEstimator == 1)
      dWeight= Cauchy::Weight(dErrorSq, dSigmaSquared);
    else
      dWeight= Huber::Weight(dErrorSq, dSigmaSquared);

    // Inlier/outlier accounting, only really works for cut-off estimators such as Tukey.
    if(dWeight < 0.5)
    {
      if(bMarkOutliers)
        TD.Point.nMEstimatorOutlierCount++;
      TD.setOutlier(true);
      continue;
    }
    else
    {
      if(bMarkOutliers)
        TD.Point.nMEstimatorInlierCount++;
      TD.setOutlier(false);
    }

    // The points are undistorded before introduce them to the estimation.    
    Vector<2> v2FoundUndist = mCamera.UnDistort(TD.v2Found);
    cv::Point2f point1;
    point1.x = v2FoundUndist[0];
    point1.y = v2FoundUndist[1];
    projPoints.push_back(point1);

    Vector<3> v3Point = TD.Point.v3WorldPos;
    cv::Point3f point2;
    point2.x = v3Point[0];
    point2.y = v3Point[1];
    point2.z = v3Point[2];
    modelPoints.push_back(point2);

  }

  cv::Mat distord = cv::Mat::zeros(5,1,CV_64F);
  cv::Mat r = cv::Mat::zeros(3,3,CV_64F);
  cv::Mat rvec = cv::Mat::zeros(3,1,CV_64F);
  cv::Mat tvec = cv::Mat::zeros(3,1,CV_64F);

  Matrix<3> R = mse3Pose.get_rotation().get_matrix();
  Vector<3> T = mse3Pose.get_translation();

  for (unsigned int i = 0; i < 3; i++)
  {
    double tmp;
    for (unsigned int j = 0; j < 3; j++)
    {
      tmp = R[i][j];
      *(r.ptr<double>(i,j)) = -tmp;
    }
    tmp = T[i];
    *(tvec.ptr<double>(0,i)) = -tmp;
  }
  cv::Rodrigues(r,rvec);

  cv::Mat & cameraMatrix = mCamera.getocvCameraMatrix();
  SE3<> pose = mse3Pose;

  if(modelPoints.size()>6 && projPoints.size()>6)
  {
    cv::solvePnP(modelPoints, projPoints, cameraMatrix, distord, rvec, tvec, bUsePrevPose);
    modelPoints.clear();
    projPoints.clear();
    distord.release();
    r.release();
    r.create(3,3,CV_64F);
    cv::Rodrigues(rvec, r);
    for(unsigned int i=0;i<3;i++)
    {
      for (unsigned int j=0;j<3;j++)
        R[i][j] = *(r.ptr<double>(i,j));
      T[i] = *(tvec.ptr<double>(0,i));
    }

    //check if we need to invert the estimation from openCV
    if((mse3Pose.get_translation()[2]>0 && T[2]<0) ||
       (mse3Pose.get_translation()[2]<0 && T[2]>0) ||
       (!bUsePrevPose && T[2]<0))
    {
      for(unsigned int i=0;i<3;i++)
      {
        for (unsigned int j=0;j<3;j++)
        {
          R[i][j] = -R[i][j];
        }
      }
      T = -1*T;
    }

    pose.get_rotation() = R;
    pose.get_translation() = T;
  }

  distord.release();
  r.release();
  rvec.release();
  tvec.release();

  return pose;
}

/**
 * This function estimates the coefficient update from the mapPoint tracking
 * information, given that the deformation basis are available.
 * 
 * @param vTD vector of pointers to NRTrackerData struct, that contains the
 *        non-rigid tracking information.
 * @param dOverrideSigma The Sigma value from which the estimator is overriden.
 * @param bMarkOutliers whether the outliers are marked or not.
 */
void EstimationUtils::CalcCoeffUpdate(std::vector<NRTrackerData *> &vTD,
                                      double dOverrideSigma,
                                      bool bMarkOutliers)
{
  if(!mpBases || !mnBasis)
    return;

  int nEstimator = 0;
  static gvar3<string> gvsEstimator("TrackerMEstimator", "Tukey", SILENT);
  if(*gvsEstimator == "Tukey")
    nEstimator = 0;
  else if(*gvsEstimator == "Cauchy")
    nEstimator = 1;
  else if(*gvsEstimator == "Huber")
    nEstimator = 2;
  else
  {
    cout << "Invalid TrackerMEstimator, choices are Tukey, Cauchy, Huber" << endl;
    nEstimator = 0;
    *gvsEstimator = "Tukey";
  }

  Vector<Dynamic,double,Reference> vDefCoefs(mpdCoeffs,mnBasis); //define the encapsulation for the coefficient set
  string strmodel = GV3::get<string>("rigid","");

  // Find the covariance-scaled reprojection error for each measurement.
  // Also, store the square of these quantities for M-Estimator sigma squared estimation.
  vector<double> vdErrorSquared, vdErrorSquared2;
  for(unsigned int f=0; f<vTD.size(); f++) // to be changed or select the indexes for the bases needed. IMPORTANT
  {
    NRTrackerData *TD = vTD[f];
    if(!TD->bFound)
      continue;
    vdErrorSquared.push_back(TD->v2Error_CovScaled * TD->v2Error_CovScaled);
  }

  //copying results to avoid being interfered by the sorting on FindSigmaSquared
  vdErrorSquared2 = vdErrorSquared; 

  // Not enough valid measurements (n>=k)? Return null update,
  // which means, not to update the deformation state.
  // There are some possible scenarios:
  // 1 - Rigid shape is available
  // 1.1 - Sequence initialized with not enough tracks to be tracked.
  //   Maintain the rigid shape coeffs, so vDefCoefs is 0, as normally is
  //   initialized
  // 1.2 - Sequence initialized with enough tracks but the current one
  //   doesn't have enough to carry on, so leave the previous state of
  //   deformation, which is the previous difference between the rigid shape
  //   and the previous one.
  // 2 - Rigid shape is not available
  // 2.1 Sequence initialized with not enough tracks to be tracked
  //   Since it is not initialized properly in any way, leave it as 0.
  // 2.2 Sequence initialized with enough tracks, but the current one
  //   doesn't have enough. Leave the previous state of deformation.
  if(vdErrorSquared.size() < mnBasis)
  {
    return;
  }

  // What is the distribution of errors?
  double dSigmaSquared;
  if(dOverrideSigma > 0)
    dSigmaSquared = dOverrideSigma;
  else
  {
    if (nEstimator == 0)
      dSigmaSquared = Tukey::FindSigmaSquared(vdErrorSquared); //quality of the estimation.
    else if(nEstimator == 1)
      dSigmaSquared = Cauchy::FindSigmaSquared(vdErrorSquared);
    else
      dSigmaSquared = Huber::FindSigmaSquared(vdErrorSquared);
  }

  WLS<Dynamic,double,SQSVD> wls(mnBasis);

  wls.add_prior(100.0); // Stabilising prior

  Vector<3> T;
  T = mse3Pose.get_translation();
  Matrix<3> R;
  R = mse3Pose.get_rotation().get_matrix();
  int a = 0;
  for(unsigned int i=0;i<vTD.size();i++) //there should be a correspondence between points and basis
  {
    NRTrackerData *TD = vTD[i];
    if(!TD->bFound)
      continue;

    double dWeight=1.0;

    if(nEstimator == 0)
      dWeight= Tukey::Weight(vdErrorSquared2[a], dSigmaSquared);
    else if(nEstimator == 1)
      dWeight= Cauchy::Weight(vdErrorSquared2[a], dSigmaSquared);
    else
      dWeight= Huber::Weight(vdErrorSquared2[a], dSigmaSquared);
     // Inlier/outlier accounting, only really works for cut-off estimators such as Tukey.

    if(dWeight < 0.5)
    {
      // Be carefull with this as, in case the map maker is running, could lead to remove the point from the map.
      if(bMarkOutliers)
        TD->Point.nMEstimatorOutlierCount++;
      a++;
      continue;
    }
    else
      if(bMarkOutliers)
        TD->Point.nMEstimatorInlierCount++;

    // The undistorded point is used to make the calculations according to the equations
    Vector<2> v2FoundUndist = mCamera.UnDistort(TD->v2Found);
    Vector<2> b_i = mCamera.GetRestTermForCoeffs(T, v2FoundUndist);
    Matrix<> a_i(2,mnBasis);
    for(unsigned int k=0; k<mnBasis; k++)
    {
      Vector<3> b = mpBases->getBasisPoint(k, a);
      Vector<2> a_ik = mCamera.GetBasisTermForCoeffs(b, v2FoundUndist, R);
      a_i.slice(0, k, 2, 1) = a_ik.as_col();
    }
    wls.add_mJ(b_i[0], a_i[0], dWeight);
    wls.add_mJ(b_i[1], a_i[1], dWeight);
    a++;
  }
  wls.compute();

  vDefCoefs = wls.get_mu();

  //take into account that if the rigid shape is loaded, the computed value for the coefficients must be taken out from those ones
  if(!strmodel.empty()){
    Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpdRigidCoeffs, mnBasis);
    vDefCoefs-=vDefCoefs_rigid;
  }

}

/**
 * This function takes the rigid shape directly from the rigid shape loaded
 * into the DBases class
 */
void EstimationUtils::rigidShapeCoeffs()
{
  if(!mpBases->isRigidFromFile())
  {
      return;
  }

  unsigned int npoints = mpBases->getNPoints();

  Matrix<> & tmp = *(mpBases->mmRigidShape);

  Matrix<> A(3*npoints, mnBasis);
  Vector<> B(3*npoints);
  A = Zeros;B = Zeros;
  for(unsigned int j=0;j<npoints;j++)
  {
    for(unsigned int i=0;i<mnBasis;i++)
    {
      A.T()[i].slice(j*3,3) = mpBases->getBasisPoint(i,j);
    }
    B.slice(j*3,3) = tmp[j];
  }
  SVD<> svd(A);
  Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpdRigidCoeffs,mnBasis);
  vDefCoefs_rigid = svd.backsub(B);
}

void EstimationUtils::rigidShapeCoeffs(Map &map)
{
  //solve the problem using least squares given all the 3D correspondences
  //Create matrix A from the already loaded bases.
  unsigned int npoints = mpBases->getNPoints();
  Vector<3> avg = Zeros;

  for(unsigned int j=0;j<npoints;j++)
  {
    avg += map.vpPoints[j]->pTData->v3Rigid;
  }

  avg /= npoints;

  for(unsigned int j=0;j<npoints;j++)
  {
    map.vpPoints[j]->pTData->v3Rigid -= avg;
  }

  Matrix<> A(3*npoints,mnBasis);
  Vector<> B(3*npoints);
  A=Zeros;
  B=Zeros;

  for(unsigned int j=0;j<npoints;j++)
  {
    for(unsigned int i=0;i<mnBasis;i++)
    {
      A.T()[i].slice(j*3,3) = mpBases->getBasisPoint(i,j)/*.as_col()*/;
    }
    B.slice(j*3,3) = map.vpPoints[j]->pTData->v3Rigid;
  }
  SVD<> svd(A);
  Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpdRigidCoeffs,mnBasis);
  vDefCoefs_rigid = svd.backsub(B);
}

/**
 * This function estimates the coefficient update from the mapPoint tracking
 * information, given that the deformation basis are available.
 *
 * @param vTD
 * @param dOverrideSigma
 * @param bMarkOutliers
 * @param base
 */
void EstimationUtils::CalcCoeffUpdate(std::vector<NRTrackerData *> &vTD,
                                      Interpolation * interpolation,
                                      double dOverrideSigma,
                                      bool bMarkOutliers)
{
  if(!interpolation || !mpBases || !mnBasis)
    return;

  //Assume that the pose is given. No initial coeffs are given
  //check what is less costy: by comparing the number of points: #base vs #detected
  interpolation->Wrap(vTD);

  int nEstimator = 0;
  static gvar3<string> gvsEstimator("TrackerMEstimator", "Tukey", SILENT);
  if(*gvsEstimator == "Tukey")
    nEstimator = 0;
  else if(*gvsEstimator == "Cauchy")
    nEstimator = 1;
  else if(*gvsEstimator == "Huber")
    nEstimator = 2;
  else
  {
    cout << "Invalid TrackerMEstimator, choices are Tukey, Cauchy, Huber" << endl;
    nEstimator = 0;
    *gvsEstimator = "Tukey";
  }

  //define the encapsulation for the coefficient set
  Vector<Dynamic,double,Reference> vDefCoefs(mpdCoeffs,mnBasis);

  // Find the covariance-scaled reprojection error for each measurement.
  // Also, store the square of these quantities for M-Estimator sigma squared estimation.
  vector<double> vdErrorSquared, vdErrorSquared2;
  for(unsigned int f=0; f<vTD.size(); f++)
  {
    NRTrackerData *TD = vTD[f];
    if(!TD->bFound ||(TD->bFound && TD->nSearchLevel>1))
      continue;
    vdErrorSquared.push_back(TD->v2Error_CovScaled * TD->v2Error_CovScaled);
  }

  vdErrorSquared2 = vdErrorSquared; //copying results to avoid later computation sorting

  // Further explanation of this condition, in the 1st implementation of this function, above
  if(vdErrorSquared.size() < mnBasis)
  {
    return;
  }

  // What is the distribution of errors?
  double dSigmaSquared;
  if(dOverrideSigma > 0)
    dSigmaSquared = dOverrideSigma;
  else
  {
    if (nEstimator == 0)
      dSigmaSquared = Tukey::FindSigmaSquared(vdErrorSquared); //quality of the estimation.
    else if(nEstimator == 1)
      dSigmaSquared = Cauchy::FindSigmaSquared(vdErrorSquared);
    else
      dSigmaSquared = Huber::FindSigmaSquared(vdErrorSquared);
  }

  WLS<Dynamic,double,SQSVD> wls(mnBasis);
  wls.add_prior(100.0);
  SE3<> se3Pose2 = mse3Pose;
  Vector<3> T;
  T = se3Pose2.get_translation();
  Matrix<3> R;
  R = se3Pose2.get_rotation().get_matrix();
  unsigned int a = 0;
  unsigned int b = 0;

  //there should be a correspondence between boints and basis
  for(unsigned int i=0;i<vTD.size();i++)
  {
    NRTrackerData *TD = vTD[i];
    if(!TD->bFound || (TD->bFound && TD->nSearchLevel>1))
      continue;

    double dWeight=1.0;

    if(nEstimator == 0)
      dWeight= Tukey::Weight(vdErrorSquared2[a], dSigmaSquared);
    else if(nEstimator == 1)
      dWeight= Cauchy::Weight(vdErrorSquared2[a], dSigmaSquared);
    else
      dWeight= Huber::Weight(vdErrorSquared2[a], dSigmaSquared);

    a++;

    // Inlier/outlier accounting, only really works for cut-off estimators such as Tukey.
    if(dWeight == 0.0)
    {
      // Be carefull with this, as explained above.
      if(bMarkOutliers)
        TD->Point.nMEstimatorOutlierCount++;
      continue;
    }
    else
      if(bMarkOutliers)
        TD->Point.nMEstimatorInlierCount++;

    //not properly initialized
    if(TD->Base==NULL)
    {
      continue;
    }

    b++;

    Vector<2> v2Point = TD->v2Found;
    //Here the undistorded point is used to make the calculations according to the equations
    Vector<2> v2FoundUndist = mCamera.UnDistort(v2Point);
    Vector<2> b_i = mCamera.GetRestTermForCoeffs(T, v2FoundUndist);
    Matrix<> a_i(2,mnBasis);
    for(unsigned int k=0; k<TD->vv3BasesPoint.size(); k++)
    {
      Vector<3> b = TD->vv3BasesPoint[k];
      Vector<2> a_ik = mCamera.GetBasisTermForCoeffs(b, v2FoundUndist, R);
      a_i.slice(0, k, 2, 1) = a_ik.as_col();
    }
    wls.add_mJ(b_i[0], a_i[0], dWeight);
    wls.add_mJ(b_i[1], a_i[1], dWeight);

    TD->updateDef3D = true;
  }

  //if there are not enough valid points for the estimation, do no compute anything else
  if(a<mnBasis)
      return;

  wls.compute();

  vDefCoefs = wls.get_mu();

  Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpdRigidCoeffs, mnBasis);

  if(vDefCoefs[0] != vDefCoefs[0])
  {
    // nan!!! not updating. setting shape to the rigid one. non-rigid to 0's
    vDefCoefs = vDefCoefs_rigid;
    return;
  }

  bool nullSolution = true;
  for(unsigned int i=0; i<mnBasis && nullSolution; i++)
    nullSolution &= (vDefCoefs[i] == 0.0);

  if(nullSolution)
  {
    vDefCoefs = vDefCoefs_rigid;
    //there should be a correspondence between boints and basis
    for(unsigned int i=0;i<vTD.size();i++)
    {
      NRTrackerData *TD = vTD[i];
      if(!TD->bFound || (TD->bFound && TD->nSearchLevel>1))
        continue;
      TD->updateDef3D = false;
    }
  }

  vDefCoefs -= vDefCoefs_rigid;

  interpolation->UnWrap(vTD);

}

/**
 * This function estimates the coefficient update from the mapPoint tracking
 * information, given that the deformation basis are available.
 *
 * @param vTD The vector of pointers to NRTrackerData, which contains Non-Rigid
 *        tracking data.
 * @param dOverrideSigma The Sigma value from which the computed one is
 *        overriden.
 * @param bMarkOutliers Whether the outliers should be marked or not.
 * @param bFirstFrame Whether the current computation is done over the first
 *        frame.
 * @param map, the current 3D map.
 * @param prevCoeff The pointer to the previous frame coefficients.
 */
void EstimationUtils::CalcCoeffUpdate2(std::vector<NRTrackerData *> &vTD,
                                       Interpolation * interpolation,
                                       double dOverrideSigma,
                                       bool bMarkOutliers,
                                       bool bFirstFrame,
                                       Map *map,
                                       double * prevCoeff)
{
  if(!interpolation || !mpBases || !mnBasis)
    return;

  // Assume that the pose is given. No initial coeffs are given
  // check what is less costy: by comparing the number of points: #base vs #detected
  interpolation->Wrap(vTD);

  bool tsmoothprior = (GV3::get<int>("timesmoothness","0")==1) & !bFirstFrame;
  double tsmoothFactor = GV3::get<double>("tsmoothfactor","1.0");

  bool shapesmoothprior = (GV3::get<int>("shapesmoothness","0")==1) & !bFirstFrame;
  double ssmoothFactor = GV3::get<double>("ssmoothfactor","1.0");

  int nEstimator = 0;
  static gvar3<string> gvsEstimator("TrackerMEstimator", "Tukey", SILENT);
  if(*gvsEstimator == "Tukey")
    nEstimator = 0;
  else if(*gvsEstimator == "Cauchy")
    nEstimator = 1;
  else if(*gvsEstimator == "Huber")
    nEstimator = 2;
  else
  {
    cout << "Invalid TrackerMEstimator, choices are Tukey, Cauchy, Huber" << endl;
    nEstimator = 0;
    *gvsEstimator = "Tukey";
  }

  Vector<Dynamic,double,Reference> vDefCoefs(mpdCoeffs,mnBasis); //define the encapsulation for the coefficient set
  // Find the covariance-scaled reprojection error for each measurement.
  // Also, store the square of these quantities for M-Estimator sigma squared estimation.
  vector<double> vdErrorSquared, vdErrorSquared2;
  for(unsigned int f=0; f<vTD.size(); f++) // to be changed or select the indexes for the bases needed. IMPORTANT
  {
    NRTrackerData *TD = vTD[f];
    if(!TD->bFound ||(TD->bFound && TD->nSearchLevel>1))
      continue;
    vdErrorSquared.push_back(TD->v2Error_CovScaled * TD->v2Error_CovScaled);
  }

  vdErrorSquared2 = vdErrorSquared; //copying results to avoid later computation

  // Further explanation of this condition, in the 1ºst implementation of this function, above
  if(vdErrorSquared.size() < mnBasis)
  {
    return;
  }

  // What is the distribution of errors?
  double dSigmaSquared;
  if(dOverrideSigma > 0)
    dSigmaSquared = dOverrideSigma; // Bit of a waste having stored the vector of square errors in this case!
  else
  {
    if (nEstimator == 0)
      dSigmaSquared = Tukey::FindSigmaSquared(vdErrorSquared); //quality of the estimation.
    else if(nEstimator == 1)
      dSigmaSquared = Cauchy::FindSigmaSquared(vdErrorSquared);
    else
      dSigmaSquared = Huber::FindSigmaSquared(vdErrorSquared);
  }

  WLS<Dynamic,double,SQSVD> wls(mnBasis);
  wls.add_prior(100.0);
  // to be continued from here to implement wls with all of these equations!!!!!!!
  //try out the different combinations
  SE3<> se3Pose2 = mse3Pose;
  Vector<3> T;
  T = se3Pose2.get_translation();
  Matrix<3> R;
  R = se3Pose2.get_rotation().get_matrix();
  unsigned int a = 0;
  unsigned int b = 0;

  for(unsigned int i=0;i<vTD.size();i++) //there should be a correspondence between boints and basis
  {
    NRTrackerData *TD = vTD[i];
    double dWeight=1.0;
    if(!TD->bFound || (TD->bFound && TD->nSearchLevel>1)){
        TD->v2Image = TD->v2Found;
        TD->v2Error_CovScaled = Zeros;
        dWeight = 0.25; // low the importance of the measurement, but not too much
    }

    if(nEstimator == 0)
      dWeight= Tukey::Weight(vdErrorSquared2[a], dSigmaSquared);
    else if(nEstimator == 1)
      dWeight= Cauchy::Weight(vdErrorSquared2[a], dSigmaSquared);
    else
      dWeight= Huber::Weight(vdErrorSquared2[a], dSigmaSquared);

    a++;

    // Inlier/outlier accounting, only really works for cut-off estimators such as Tukey.
    if(dWeight < 0.4)
    {
      // Caution, if the mapmaker thread is active it could yield in
      // map degradation as it would erase the point.
      if(bMarkOutliers)
        TD->Point.nMEstimatorOutlierCount++;
      continue;
    }
    else
      if(bMarkOutliers)
        TD->Point.nMEstimatorInlierCount++;
    //not properly initialized
    if(TD->Base==NULL) 
    {
      continue;
    }

    b++;

    Vector<2> v2Point = TD->v2Found;
    // Here the undistorded point is used to make the calculations according to the equations
    Vector<2> v2FoundUndist = mCamera.UnDistort(v2Point);
    Vector<2> b_i = mCamera.GetRestTermForCoeffs(T, v2FoundUndist, TD->v3Rigid, R);
    Matrix<> a_i(2,mnBasis);
    Matrix<> tprior(3,mnBasis);
    for(unsigned int k=0; k<TD->vv3BasesPoint.size(); k++)
    {
      Vector<3> b = TD->vv3BasesPoint[k];
      Vector<2> a_ik = mCamera.GetBasisTermForCoeffs(b, v2FoundUndist, R);
      a_i.slice(0, k, 2, 1) = a_ik.as_col();
      if (tsmoothprior)
      {
        tprior.T()[k] = b;
      }
    }
    wls.add_mJ(b_i[0], a_i[0], dWeight);
    wls.add_mJ(b_i[1], a_i[1], dWeight);

    if (tsmoothprior){
      wls.add_mJ(TD->v3MapPrev[0], tprior[0],tsmoothFactor*dWeight);
      wls.add_mJ(TD->v3MapPrev[1], tprior[1],tsmoothFactor*dWeight);
      wls.add_mJ(TD->v3MapPrev[2], tprior[2],tsmoothFactor*dWeight);
    }
    TD->updateDef3D = true;
  }
  
  if (shapesmoothprior){
    // let's assume the indices of the map is the same as the ones on the bases
    // and the same on the m3Shape on the interpolation object
    Matrix<> sprior(3,mnBasis);
    for (unsigned int j = 0; j < interpolation->m3Shape.size(); j++){
      vector<unsigned int> &vidx = interpolation->getNeighboursFromPointIndex(j);
      Vector<3,unsigned int> selidx = interpolation->getSelectedNeighboursIdx(j);
      if (vidx.size() < 3 || (selidx[0]==0 && selidx[1]==0 && selidx[2]==0))
        continue;
      Vector<3> ct = interpolation->getConstraintsNeighbours(j);
      Vector<3> c = interpolation->m3Shape[j];
      for (unsigned int k = 0; k < mnBasis; k++){
        Vector<3> b1 = mpBases->getBasisPoint(k,vidx[selidx[0]]);
        Vector<3> b2 = mpBases->getBasisPoint(k,vidx[selidx[1]]);
        Vector<3> b3 = mpBases->getBasisPoint(k,vidx[selidx[2]]);
        sprior.T()[k] = (ct[0]*b1+ct[1]*b2+ct[2]*b3);
      }
      wls.add_mJ(c[0], sprior[0], ssmoothFactor);
      wls.add_mJ(c[1], sprior[1], ssmoothFactor);
      wls.add_mJ(c[2], sprior[2], ssmoothFactor);    
    }
  }
  //if there are not enough valid points for the estimation, do not compute anything else
  if(a<mnBasis)
      return;

  wls.compute();

  vDefCoefs = wls.get_mu();

  Vector<Dynamic,double,Reference> vDefCoefs_rigid(mpdRigidCoeffs, mnBasis);

  if(vDefCoefs[0] != vDefCoefs[0])
  {
    // Nan!!! not updating. setting the shape to the rigid one. non-rigid to 0's
    vDefCoefs = vDefCoefs_rigid;
    return;
  }

  bool nullSolution = true;
  for(unsigned int i=0; i<mnBasis && nullSolution; i++)
    nullSolution &= (vDefCoefs[i] == 0.0);

  if(nullSolution)
  {
    vDefCoefs = vDefCoefs_rigid;
    //there should be a correspondence between boints and basis
    for(unsigned int i=0;i<vTD.size();i++)
    {
      NRTrackerData *TD = vTD[i];
      if(!TD->bFound || (TD->bFound && TD->nSearchLevel>1))
        continue;
      TD->updateDef3D = false;
    }
  }

  vDefCoefs -= vDefCoefs_rigid;
  interpolation->UnWrap(vTD);
}
