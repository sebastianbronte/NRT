/*
 * VisualNRTracker.h
 * Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
 * Date 13/01/2013
 * 
 * This class performs specific tracking for visual input for Non-Rigid tracks.
 */

#ifndef __VISUALNRTRACKER__
#define __VISUALNRTRACKER__

#include "Tracker.h"
#include "EstimationUtils.h"
#include "DBases.h"
#include <opencv/cv.h>
#include "DescriptorMatcherOCV.h"

class Interpolation;

class VisualNRTracker : public Tracker
{
  public:
    VisualNRTracker(CVD::ImageRef irVideoSize, const ATANCamera &c, Map &m,
                    DBases &base);
    virtual ~VisualNRTracker();

    void Reset();

    virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame, bool bDraw);

    virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame,
                            std::vector<double> &coeff,
                            bool bDraw);

    virtual void TrackFrame(std::vector<Vector<2> > &tracks,
                          std::vector<bool> &visibility,
                          std::vector<double> &coeffs,
                          bool bDraw=true){};

    virtual bool AreTracksLoadedFromDisk(){return false;}

    virtual unsigned int getNBasis(){return mBase.getNBasis();}
    virtual unsigned int getNBasisPoints(){return mBase.getNPoints();}
    virtual const double * getCoefs(){return mpDefCoefs;}

    virtual std::vector<NRTrackerData*>& vIterationSet() {return vIterationSetNR;}

    virtual Interpolation* getInterpolation() {return mpInterpolation;}

  protected:
    virtual void TrackMap();

    virtual void TrackForInitialMap();

    int SearchForPoints(std::vector<NRTrackerData*> &vTD,
                        int nRange, int nSubPixIts);

    std::vector<NRTrackerData*> vIterationSetNR;
    //used to perform trail tracking on the map
    std::map<NRTrackerData *, Trail> mMapTrail;

    ///Deformation stuff

    //Deformation bases
    DBases & mBase;
    //dynamic coefficients
    double * mpDefCoefs;
    // In case a rigid shape can be expressed as a combination of the
    // basis shapes
    double * mpDefCoefs_rigid;
    // The last known working set of coefficients
    double * mpDefCoefs_default;

    //Descriptor matcher based on descriptor
    DescriptorMatcherOCV mdMatcher;

    //Interpolation needed to estimate the baricentric coordinates
    //of the detected points
    Interpolation *mpInterpolation;

    // Overload of the function on the base class so as to call it with
    // parameters enough for priors
    Vector<6> CalcPoseUpdate(std::vector<NRTrackerData*> &vTD,
                         double dOverrideSigma = 0.0,
                         bool bMarkOutliers = false);

    // Coefficient estimation function. Calls the corresponding EstimationUtils
    // function
    void CalcCoeffUpdate(std::vector<NRTrackerData*> &vTD,
                         double dOverrideSigma = 0.0,
                         bool bMarkOutliers = false);

    void setUpdate3DFlags(std::vector<NRTrackerData*> &vTD, bool value=true);
    void bootstrapDeformation();

    void alignMapFromDisk();

    //void expandMap();

    // Priors: maintain some data from the previous shape to apply time smoothness
    SE3<> mse3CamFromWorldPrev;
    // the best last frame coefficients
    double * mpDefCoefsPrev;
};

#endif
