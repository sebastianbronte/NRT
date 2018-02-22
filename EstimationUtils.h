/**
 * EstimationUtils.h
 * Author: Sebastián Bronte <sebastian.bronte@depeca.uah.es>
 * Date: 13/01/2013
 * 
 * Class to implement the different estimation functions, making them 
 * independent from the Tracking, and eventually Mapping, process
 * 
 */

#ifndef __ESTIMATIONUTILS__
#define __ESTIMATIONUTILS__

#include "ATANCamera.h"
#include "Relocaliser.h"

#include <vector>

class ATANCamera;
class NRTrackerData;
class DBases;
class Map;
class Interpolation;

class EstimationUtils
{
  public:

    EstimationUtils(ATANCamera &camera, SE3<> & pose);

    void setCoefficients(double * pCoeffs, double * pRigidCoeffs,
                         unsigned int nBasis);
    void setBases(DBases * bases);
    void setPose(SE3<> & se3CamFromPose);

    // This version computes the pose update based on jacobian computation
    // and least squares
    Vector<6> CalcPoseUpdate(std::vector<NRTrackerData *> &vTD,
                             double dOverrideSigma = 0.0,
                             bool bMarkOutliers = false,
                             bool bFirstFrame=false);

    void CalcCoeffUpdate(std::vector<NRTrackerData *> &vTD,
                         double dOverrideSigma = 0.0,
                         bool bMarkOutliers = false);

    void CalcCoeffUpdate(std::vector<NRTrackerData *> &vTD,
                         Interpolation * interpolation,
                         double dOverrideSigma = 0.0,
                         bool bMarkOutliers = false);

    void CalcCoeffUpdate2(std::vector<NRTrackerData *> &vTD,
                          Interpolation * interpolation,
                          double dOverrideSigma = 0.0,
                          bool bMarkOutliers = false,
                          bool bfirstframe=false,
                          Map *map=NULL,
                          double * prevCoeff=NULL);

    SE3<> CalcPoseEPnP(std::vector<NRTrackerData *> &vTD,
                       double dOverrideSigma,
                       bool bMarkOutliers,
                       bool bUsePrevPose);

    void rigidShapeCoeffs(Map &map);
    void rigidShapeCoeffs();
  protected:
    ATANCamera &mCamera;
    SE3<> &mse3Pose;
    double * mpdCoeffs;
    double * mpdRigidCoeffs;
    unsigned int mnBasis;

    DBases *mpBases;
};

#endif