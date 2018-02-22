//-*- C++ -*-
// Copyright 2008 Isis Innovation Limited
// 
// This header declares the Tracker class.
// The Tracker is one of main components of the system,
// and is responsible for determining the pose of a camera
// from a video feed. It uses the Map to track, and communicates 
// with the MapMaker (which runs in a different thread)
// to help construct this map.
//
// Initially there is no map, so the Tracker also has a mode to 
// do simple patch tracking across a stereo pair. This is handled 
// by the TrackForInitialMap() method and associated sub-methods. 
// Once there is a map, TrackMap() is used.
//
// Externally, the tracker should be used by calling TrackFrame()
// with every new input video frame. This then calls either 
// TrackForInitialMap() or TrackMap() as appropriate.
//

#ifndef __TRACKER_H
#define __TRACKER_H

#include "EstimationUtils.h"
#include "KeyFrame.h"
#include "ATANCamera.h"
#include "MiniPatch.h"
#include "Relocaliser.h"

#include <sstream>
#include <vector>
#include <list>

#define DEBUG false

#define NUM_LOST_FRAMES 3

class Interpolation;

class NRTrackerData;

struct Trail    // This struct is used for initial correspondences of the first stereo pair.
{
  MiniPatch mPatch;
  CVD::ImageRef irCurrentPos;
  CVD::ImageRef irInitialPos;
  int nSearchLevel;
};

class AbstractTracker
{
public:
  virtual ~AbstractTracker(){}
  virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame, bool bDraw) = 0;
  virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame,
                          std::vector<double> &coeff,
                          bool bDraw) = 0;
  virtual void TrackFrame(std::vector<Vector<2> > &tracks,
                          std::vector<bool> &visibility,
                          std::vector<double> &coeffs,
                          bool bDraw=true) = 0;

  virtual int GetFrameNumber() = 0;

  virtual std::string GetMessageForUser() = 0;

  virtual SE3<> GetCurrentPose() = 0;

  virtual bool AreCoeffsLoadedFromDisk() = 0;
  virtual bool AreTracksLoadedFromDisk() = 0;
  virtual std::vector<NRTrackerData*>& getvIterationSet() = 0;
  virtual int getTrackingQuality() = 0;
  virtual unsigned int getNBasis() = 0;
  virtual const double * getCoefs() = 0;
  virtual Interpolation* getInterpolation() = 0;
  virtual KeyFrame & getCurrentKF() = 0;
protected:
  virtual void TrackMap() = 0;
  virtual void TrackForInitialMap() = 0;
  virtual void TrailTracking_Start() = 0;
  virtual int TrailTracking_Advance() = 0;
};

class Tracker : public AbstractTracker
{
public:
  Tracker(CVD::ImageRef irVideoSize, const ATANCamera &c, Map &m, int levels = 4);
  virtual ~Tracker();

  // TrackFrame is the main working part of the tracker: call this every frame.
  virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame, bool bDraw);
  // This version is called when the pre-defined coefficients or pose are used.
  virtual void TrackFrame(CVD::Image<CVD::byte> &imFrame,
                          std::vector<double> &coeff,
                          bool bDraw){};

  // Alternative to tracking if the tracking algorithm is point-wise
  virtual void TrackFrame(std::vector<Vector<2> > &tracks,
                          std::vector<bool> &visibility,
                          std::vector<double> &coeffs,
                          bool bDraw=true){}

  virtual inline SE3<> GetCurrentPose() { return mse3CamFromWorld;}

  // Gets messages to be printed on-screen for the user.
  virtual std::string GetMessageForUser();

  virtual std::vector<NRTrackerData*>& getvIterationSet() { return vIterationSet; }

  virtual int getTrackingQuality() { return mTrackingQuality; }

  virtual inline int GetFrameNumber() {return mnFrame;}

  virtual bool AreTracksLoadedFromDisk(){ return false;}
  virtual bool AreCoeffsLoadedFromDisk(){ return areCoeffsLoadedFromDisk; }

  virtual unsigned int getNBasis() { return 0; }
  virtual const double * getCoefs() { return NULL; }

  virtual Interpolation* getInterpolation() {return NULL;}

  virtual KeyFrame & getCurrentKF() { return mCurrentKF; }

protected:
  KeyFrame mCurrentKF;            // The current working frame as a keyframe struct
  // The major components to which the tracker needs access:
  Map &mMap;                      // The map, consisting of points and keyframes
  ATANCamera mCamera;             // Projection model
  Relocaliser mRelocaliser;       // Relocalisation module

  CVD::ImageRef mirSize;          // Image size of whole image

  void Reset();                   // Restart from scratch. Also tells the mapmaker to reset itself.
  void RenderGrid();              // Draws the reference grid

  enum {BAD, DODGY, GOOD} mTrackingQuality;

  // The following members are used for initial map tracking (to get the first stereo pair and correspondences):
  virtual void TrackForInitialMap();      // This is called by TrackFrame if there is not a map yet.
  enum {TRAIL_TRACKING_NOT_STARTED, 
        TRAIL_TRACKING_STARTED, 
        TRAIL_TRACKING_COMPLETE} mnInitialStage;  // How far are we towards making the initial map?
  virtual void TrailTracking_Start();     // First frame of initial trail tracking. Called by TrackForInitialMap.
  virtual int  TrailTracking_Advance();   // Steady-state of initial trail tracking. Called by TrackForInitialMap.
  std::list<Trail> mlTrails;      // Used by trail tracking
  KeyFrame mFirstKF;              // First of the stereo pair
  KeyFrame mPreviousFrameKF;      // Used by trail tracking to check married matches
  std::vector<NRTrackerData*> vIterationSet;

  // Methods for tracking the map once it has been made:
  virtual void TrackMap();                // Called by TrackFrame if there is a map.
  void AssessTrackingQuality();   // Heuristics to choose between good, poor, bad.
  void ApplyMotionModel();        // Decaying velocity motion model applied prior to TrackMap
  void UpdateMotionModel();       // Motion model is updated after TrackMap
  int SearchForPoints(std::vector<NRTrackerData*> &vTD, 
                      int nRange, 
                      int nFineIts);  // Finds points in the image
  Vector<6> CalcPoseUpdate(std::vector<NRTrackerData*> &vTD,
                           double dOverrideSigma = 0.0, 
                           bool bMarkOutliers = false); // Updates pose from found points.
  SE3<> mse3CamFromWorld;           // Camera pose: this is what the tracker updates every frame.
  SE3<> mse3StartPos;               // What the camera pose was at the start of the frame.
  Vector<6> mv6CameraVelocity;    // Motion model
  double mdVelocityMagnitude;     // Used to decide on coarse tracking 
  double mdMSDScaledVelocityMagnitude; // Velocity magnitude scaled by relative scene depth.
  bool mbDidCoarse;               // Did tracking use the coarse tracking stage?

  bool mbDraw;                    // Should the tracker draw anything to OpenGL?

  // Interface with map maker:
  int mnFrame;                    // Frames processed since last reset
  int mnLastKeyFrameDropped;      // Counter of last keyframe inserted.
  //void AddNewKeyFrame();          // Gives the current frame to the mapmaker to use as a keyframe

  // Tracking quality control:
  int manMeasAttempted[LEVELS];
  int manMeasFound[LEVELS];
  int mnLostFrames;

  // Relocalisation functions:
  bool AttemptRecovery();         // Called by TrackFrame if tracking is lost.
  bool mbJustRecoveredSoUseCoarse;// Always use coarse tracking after recovery!

  // Frame-to-frame motion init:
  SmallBlurryImage *mpSBILastFrame;
  SmallBlurryImage *mpSBIThisFrame;
  void CalcSBIRotation();
  Vector<6> mv6SBIRot;
  bool mbUseSBIInit;

  // User interaction for initial tracking:
  bool mbUserPressedSpacebar;
  std::ostringstream mMessageForUser;
  //user-defined initial frames
  int iniframe1;
  int iniframe2;

  bool areCoeffsLoadedFromDisk;

  // GUI interface:
  void GUICommandHandler(std::string sCommand, std::string sParams);
  static void GUICommandCallBack(void* ptr, std::string sCommand, std::string sParams);
  struct Command {std::string sCommand; std::string sParams; };
  std::vector<Command> mvQueuedCommands;

  //error estimation to decide whether update the coefficients or not
  double mdRMSError;
  double computeProjError(std::vector<NRTrackerData*> &vTD);
  double computeProjError_RMS(std::vector<NRTrackerData*> &vTD);
  double computeBestKFProjError_RMS(std::vector<NRTrackerData*> &vTD);

  EstimationUtils mStateEstimator;

  bool mbMapGrowing;
};

#endif
