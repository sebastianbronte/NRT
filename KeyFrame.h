// -*- c++ -*-
// Copyright 2008 Isis Innovation Limited

//
// This header declares the data structures to do with keyframes:
// structs KeyFrame, Level, Measurement, Candidate.
// 
// A KeyFrame contains an image pyramid stored as array of Level;
// A KeyFrame also has associated map-point mesurements stored as a vector of Measurment;
// Each individual Level contains an image, corner points, and special corner points
// which are promoted to Candidate status (the mapmaker tries to make new map points from those.)
//
// KeyFrames are stored in the Map class and manipulated by the MapMaker.
// However, the tracker also stores its current frame as a half-populated
// KeyFrame struct.


#ifndef __KEYFRAME_H
#define __KEYFRAME_H
#include "SmallBlurryImage.h"
#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <cvd/image.h>
#include <cvd/byte.h>
#include <vector>
#include <set>
#include <map>
#include <gvars3/instances.h>

using namespace TooN;

class NRMapPoint;
class DescriptorMatcherOCV;

// Regular visual tracking. This indicates the amount of levels of
// the pyramid based tracking.
#define LEVELS 4
// Uncomment to perform model based tracking for the text-based tracker, apply
// the algorithm in the first pyramid level only, and comment the previous definition.
//#define LEVELS 1

// Candidate: a feature in an image which could be made into a map point
struct Candidate
{
  CVD::ImageRef irLevelPos;
  Vector<2> v2RootPos;
  double dSTScore;
};

// Measurement: A 2D image measurement of a map point. Each keyframe stores a bunch of these.
struct Measurement
{
  int nLevel;   // Which image level?
  bool bSubPix; // Has this measurement been refined to sub-pixel level?
  Vector<2> v2RootPos;  // Position of the measurement, REFERED TO PYRAMID LEVEL ZERO
  enum {SRC_TRACKER, SRC_REFIND, SRC_ROOT, SRC_TRAIL, SRC_EPIPOLAR} Source; // Where has this measurement come from?
};

// Each keyframe is made of LEVELS pyramid levels, stored in struct Level.
// This contains image data and corner points.
struct Level
{
  inline Level()
  {
    bImplaneCornersCached = false;
    vCorners.clear();
    vCornerRowLUT.clear();
    vCandidates.clear();
    vImplaneCorners.clear();
  }
  
  CVD::Image<CVD::byte> im;                // The pyramid level pixels
  std::vector<CVD::ImageRef> vCorners;     // All FAST corners on this level
  std::vector<int> vCornerRowLUT;          // Row-index into the FAST corners, speeds up access
  std::vector<CVD::ImageRef> vMaxCorners;  // The maximal FAST corners
  Level& operator=(const Level &rhs);
  
  std::vector<Candidate> vCandidates;      // Potential locations of new map points
  
  bool bImplaneCornersCached;              // Also keep image-plane (z=1) positions of FAST corners to speed up epipolar search
  std::vector<Vector<2> > vImplaneCorners; // Corner points un-projected into z=1-plane coordinates
};

// The actual KeyFrame struct. The map contains of a bunch of these. However, the tracker uses this
// struct as well: every incoming frame is turned into a keyframe before tracking; most of these 
// are then simply discarded, but sometimes they're then just added to the map.
struct KeyFrame
{
  // Set the maximum number of levels to be used for this KF (4 maximum)
  inline KeyFrame(int maxlevels = LEVELS):pSBI(NULL),MAXLEVELS(maxlevels) {}
  inline ~KeyFrame() {}
  // The coordinate frame of this key-frame as a Camera-From-World transformation
  SE3<> se3CfromW;
  // Is the coordinate frame of this keyframe fixed? (only true for first KF!)
  bool bFixed;
  // Images, corners, etc live in this array of pyramid levels
  Level aLevels[LEVELS];
  // All the measurements associated with the keyframe
  std::map<NRMapPoint*, Measurement> mMeasurements;

  // Takes an image and computes pyramid levels to fill the keyframe
  // data structures with what is needed by the tracker
  void MakeKeyFrame_Lite(CVD::BasicImage<CVD::byte> &im);

  // Takes an image and the features from an OpenCV descriptor to fill the
  // keyframe data structures with what is needed by the tracker.
  void MakeKeyFrame_Lite(const CVD::BasicImage<CVD::byte> &im, DescriptorMatcherOCV &dMatcher);   


  void MakeKeyFrame_Rest();              // ... while this calculates the rest of the data which the mapmaker needs.
  void MakeKeyFrame_Rest2();

  double dSceneDepthMean;      // Hacky hueristics to improve epipolar search.
  double dSceneDepthSigma;

  SmallBlurryImage *pSBI; // The relocaliser uses this

  int MAXLEVELS;
};

typedef std::map<NRMapPoint*, Measurement>::iterator meas_it;  // For convenience, and to work around an emacs paren-matching bug

#endif
