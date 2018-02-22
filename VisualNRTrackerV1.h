/*
 * VisualNRTracker.h
 * Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
 * Date 13/01/2013
 * 
 * This class performs specific tracking for visual input for Non-Rigid tracks.
 */

#ifndef __VISUALNRTRACKERV1__
#define __VISUALNRTRACKERV1__

#include "Tracker.h"
#include "EstimationUtils.h"
#include "DBases.h"

class Interpolation;

class VisualNRTrackerV1 : public Tracker
{
  public:
    VisualNRTrackerV1(CVD::ImageRef irVideoSize,
                      const ATANCamera &c,
                      Map &m,
                      DBases &base);
    virtual ~VisualNRTrackerV1();

    void Reset();

    virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame, bool bDraw);
    virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame,
                            std::vector<double> &coeff,
                            bool bDraw);
    //only for debug purposes
    virtual void TrackFrame(std::vector<Vector<2> > &tracks,
                            std::vector<bool> &visibility,
                            std::vector<double> &coeffs,
                            bool bDraw=true){};

    virtual bool AreTracksLoadedFromDisk(){return false;}

    virtual unsigned int getNBasis(){return mBase.getNBasis();}
    virtual unsigned int getNBasisPoints(){return mBase.getNPoints();}
    virtual const double * getCoefs(){return mpDefCoefs;}

    virtual     std::vector<NRTrackerData*>& getvIterationSet(){return vIterationSetNR;}

    virtual Interpolation* getInterpolation() {return mpInterpolation;}

  protected:
    virtual void TrackMap();

    virtual void TrackForInitialMap();
    int SearchForPoints(std::vector<NRTrackerData*> &vTD,
                        int nRange, int nSubPixIts);
    //Helper function to perform just trail tracking on the MapPoints
    bool SecondChanceSearch(NRTrackerData* TD, int nRange);

    std::vector<NRTrackerData*> vIterationSetNR;
    //used to perform trail tracking on the map
    std::map<NRTrackerData *, Trail> mMapTrail;

    //Deformation stuff

    //Deformation bases
    DBases & mBase;
    //dynamic coefficients
    double * mpDefCoefs;
    //a rigid shape could be expressed as a combination of the available set of basis.
    double * mpDefCoefs_rigid;
    // the last known working set of coefficients
    double * mpDefCoefs_default;

    //Interpolation needed to estimate the baricentric coordinates
    //of the detected points
    Interpolation *mpInterpolation;

    //Overload of the function on the base class so as to call it with parameters
    //enough for priors
    Vector<6> CalcPoseUpdate(std::vector<NRTrackerData*> &vTD,
                         double dOverrideSigma = 0.0,
                         bool bMarkOutliers = false);

    //Coefficient estimation function -> it will call the EstimationUtils 
    //function especiallized in this.
    void CalcCoeffUpdate(std::vector<NRTrackerData*> &vTD,
                         double dOverrideSigma = 0.0,
                         bool bMarkOutliers = false);

    void setUpdate3DFlags(std::vector<NRTrackerData*> &vTD, bool value=true);
    void bootstrapDeformation();

    //some stuff if the map has been loaded from disk. Temporally leave it apart. Some kind of alignment will be given
    //inherited in a way from NRTracker class
    void alignMapFromDisk();

    //priors: maintain some data from the previous shape to apply time smoothness
    SE3<> mse3CamFromWorldPrev; 
};

#endif