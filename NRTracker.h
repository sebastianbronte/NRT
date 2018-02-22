// -*- c++ -*-
// Author: Sebastian Bronte Palacios. sebastian.bronte@depeca.uah.es
// Date: 3 Feb 2011
// NRTracker.h
//
// Defines the NRTracker class and aditional structures needed: DBase and NRTrackerData
//
// This class handles some algorithms to perform NR-SfM on PTAM. In this file is
// only included the modified version of the original tracker, to handle
// deforming objects from synthetic data.
//

#ifndef _NR_TRACKER_
#define _NR_TRACKER_

#include "Tracker.h"
#include "ATANCamera.h"
#include "MiniPatch.h"
#include "Relocaliser.h"
#include "DBases.h"

#include <sstream>
#include <vector>
#include <list>


class NRTrackerData;
class EstimationUtils;
class Interpolation;

/**
 * This class is written in order to substitute the classical tracking process
 * implemented in PTAM. Some functions will be overridden to adapt then to the
 * new framework.
 */

class NRTracker: public Tracker
{
public:

    NRTracker(CVD::ImageRef &irVideoSize, const ATANCamera &c, Map &m, DBases &base);
    virtual ~NRTracker();

    void Reset();

    virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame, bool bDraw);
    virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame,
                            std::vector<double> &coeff,
                            bool bDraw){}
    virtual void TrackFrame(std::vector<Vector<2> > &tracks, std::vector<bool> &visibility, std::vector<double> &coeffs, bool bDraw=true);

    virtual std::vector<NRTrackerData*>& getvIterationSet() { return vIterationSetNR; }

    virtual bool AreTracksLoadedFromDisk(){ return areTracksLoadedFromDisk; }

    inline unsigned int getNBasis(){return base.getNBasis();}
    inline unsigned int getNBasisPoints(){return base.getNPoints();}
    inline const double * getCoefs(){return mpDefCoefs;}
    virtual Interpolation * getInterpolation() {return NULL;}

protected:
    DBases &base;
    double *mpDefCoefs;
    double *mpDefCoefs_rigid;

    //The map has been loaded from disk. Some alignment will be needed.
    void alignMapFromDisk();

    //rewritten functions
    //helper function, just in case we want to perform the original or the NRSfM method.
    void TrackFrameFromGivenTracks(); 
    virtual void TrackMap();

    //deformation coefficient functions
    //update the coefficients from the current state of the program
    void CalcCoeffUpdate(std::vector<NRTrackerData*> &vTD);
    //estimate the rigid coefficients given a good pose and the basis shapes
    void rigidShapeCoeffs();

    void setUpdate3DFlags(std::vector<NRTrackerData*> &vTD);

    /// Relocalisation functions:
    // Called by TrackFrame if tracking is lost.
    bool AttemptRecovery();
    SE3<> mse3CamFromWorld_default;
    double *mpDefCoefs_default;

    std::vector<NRTrackerData*> vIterationSetNR;

    //where the tracks from file will live
    std::vector<Vector<2> > *mpTracks;

    bool areTracksLoadedFromDisk;
    bool applyMotionModel;

    bool normalizedModel;
    void normalizeCoefficients();
    void normalizeRigidCoefficients();
    // Decaying velocity deformation model applied prior to NRTrackMap
    void ApplyDeformationModelSize();
    // Deformation model is updated after NRTrackMap
    void UpdateDeformationModelSize();
    double rigidShapeEnergy;
    double currentShapeEnergy;
    double defModelUpdate;
};

#endif
