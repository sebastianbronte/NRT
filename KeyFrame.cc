// Copyright 2008 Isis Innovation Limited
#include "KeyFrame.h"
#include "ShiTomasi.h"
#include "SmallBlurryImage.h"
#include "DescriptorMatcherOCV.h"
#include "OpenCVUtilities.h"
#include <cvd/vision.h>
#include <cvd/fast_corner.h>

#include <opencv/cv.h>

using namespace CVD;
using namespace std;
using namespace GVars3;

/**
 * Perpares a Keyframe from an image. Generates pyramid levels, does FAST detection, etc.
 * Does not fully populate the keyframe struct, but only does the bits needed for the tracker;
 * e.g. does not perform FAST nonmax suppression. Things like that which are needed by the
 * mapmaker but not the tracker go in MakeKeyFrame_Rest();
 * @param im image to make keyframe from
 */
void KeyFrame::MakeKeyFrame_Lite(BasicImage<CVD::byte> &im)
{
  // First, copy out the image data to the pyramid's zero level.
  aLevels[0].im.resize(im.size());
  copy(im, aLevels[0].im);

  // Then, for each level...
  for(int i=0; i<MAXLEVELS; i++)
  {
    Level &lev = aLevels[i];
    if(i!=0)
    {  // .. make a half-size image from the previous level..
      lev.im.resize(aLevels[i-1].im.size() / 2);
      halfSample(aLevels[i-1].im, lev.im);
    }

    // .. and detect and store FAST corner points.
    // I use a different threshold on each level; this is a bit of a hack
    // whose aim is to balance the different levels' relative feature densities.
    lev.vCorners.clear();
    lev.vCandidates.clear();
    lev.vMaxCorners.clear();
    if(i == 0)
      fast_corner_detect_10(lev.im, lev.vCorners, 10);
    if(i == 1)
      fast_corner_detect_10(lev.im, lev.vCorners, 15);
    if(i == 2)
      fast_corner_detect_10(lev.im, lev.vCorners, 15);
    if(i == 3)
      fast_corner_detect_10(lev.im, lev.vCorners, 10);

    // Generate row look-up-table for the FAST corner points: this speeds up 
    // finding close-by corner points later on.
    unsigned int v=0;
    lev.vCornerRowLUT.clear();
    for(int y=0; y<lev.im.size().y; y++)
    {
      while( (v < lev.vCorners.size()) && (y > lev.vCorners[v].y) )
        v++;
      lev.vCornerRowLUT.push_back(v);
    }
  }
}

void KeyFrame::MakeKeyFrame_Lite(const BasicImage<CVD::byte> &im, DescriptorMatcherOCV &dMatcher)
{
  // First, copy out the image data to the pyramid's zero level.
  ImageRef s = im.size();
  aLevels[0].im.resize(s);
  copy(im, aLevels[0].im);
  cv::Mat ocvImage = cv::Mat(s.y, s.x, CV_8UC1);

  std::vector<cv::KeyPoint> kpts;

  OpenCVUtilities::conversionCVD2OpenCV_NB(im, ocvImage);
  dMatcher.extractFeaturesAndDescriptorsFromImage(ocvImage); //the input image is set, the features
  kpts = dMatcher.getPoints();

  for (int l=0;l<MAXLEVELS;l++)
  {
    Level &lev = aLevels[l];
    if (l!=0)
    {
      lev.im.resize(aLevels[l-1].im.size()/2);
      halfSample(aLevels[l-1].im, lev.im);
    }
    lev.vCorners.clear();
    lev.vCandidates.clear();
    lev.vMaxCorners.clear();
  }

  //now there the levels are undefined, as they will be defined by its octave
  //so the following loop compiles the data for each point
  for(unsigned int i=0;i<kpts.size();i++)
  {
    cv::KeyPoint &kpt = kpts[i];
    int l =  kpt.octave;
    if (l >= MAXLEVELS)
      l = MAXLEVELS-1; //In fact the number of levels is just one, but ... just in case.
    Level &lev = aLevels[l];

    ImageRef pt; pt.x = kpt.pt.x;pt.y = kpt.pt.y;

    lev.vCorners.push_back(pt);

  }

  //Is the LUT needed anymore? Not for tracking.
  //In case Bundle Adjustment starts it is necessary. Not in this release.
  /*for (unsigned int l = 0; l < MAXLEVELS; l++){
    unsigned int v=0;
    Level &lev = aLevels[l];
    lev.vCornerRowLUT.clear();
    for(int y=0; y<lev.im.size().y; y++)
    {
      while( (v < lev.vCorners.size()) && (y > lev.vCorners[v].y) )
        v++;
      lev.vCornerRowLUT.push_back(v);
    }
  }*/
}

/**
 * Fills the rest of the keyframe structure needed by the mapmaker:
 * FAST nonmax suppression, generation of the list of candidates for further map points,
 * creation of the relocaliser's SmallBlurryImage.
 */
void KeyFrame::MakeKeyFrame_Rest()
{
  static gvar3<double> gvdCandidateMinSTScore("MapMaker.CandidateMinShiTomasiScore", 70, SILENT);
  // For each level...
  for(int l=0; l<MAXLEVELS; l++)
  {
    Level &lev = aLevels[l];
    // .. find those FAST corners which are maximal..
    fast_nonmax(lev.im, lev.vCorners, 10, lev.vMaxCorners);
    // .. and then calculate the Shi-Tomasi scores of those, and keep the ones with
    // a suitably high score as Candidates, i.e. points which the mapmaker will attempt
    // to make new map points out of.
    for(vector<ImageRef>::iterator i=lev.vMaxCorners.begin(); i!=lev.vMaxCorners.end(); i++)
    {
      if(!lev.im.in_image_with_border(*i, 10))
        continue;
      double dSTScore = FindShiTomasiScoreAtPoint(lev.im, 3, *i);
      if(dSTScore > *gvdCandidateMinSTScore)
      {
        Candidate c;
        c.irLevelPos = *i;
        c.dSTScore = dSTScore;
        lev.vCandidates.push_back(c);
      }
    }
  }

  // Also, make a SmallBlurryImage of the keyframe: The relocaliser uses these.
  pSBI = new SmallBlurryImage(*this);  
  // Relocaliser also wants the jacobians..
  pSBI->MakeJacs();
}

void KeyFrame::MakeKeyFrame_Rest2()
{
  static gvar3<double> gvdCandidateMinSTScore("MapMaker.CandidateMinShiTomasiScore", 70, SILENT);

  for(int l=0; l<MAXLEVELS; l++)
  {
    Level &lev = aLevels[l];
    for(vector<ImageRef>::iterator i=lev.vMaxCorners.begin(); i!=lev.vMaxCorners.end(); i++)
    {
      if(!lev.im.in_image_with_border(*i, 10))
        continue;
      double dSTScore = 100;
      if(dSTScore > *gvdCandidateMinSTScore)
      {
        Candidate c;
        c.irLevelPos = *i;
        c.dSTScore = dSTScore;
        lev.vCandidates.push_back(c);
      }
    }
  }

  // Also, make a SmallBlurryImage of the keyframe: The relocaliser uses these.
  //if (pSBI) {delete pSBI; pSBI=NULL;}
  pSBI = new SmallBlurryImage(*this);  
  // Relocaliser also wants the jacobians..
  pSBI->MakeJacs();
}


/**
 * Level needs its own operator= to override CVD's reference-counting behaviour.
 * @param rhs the level to copy
 * @return the copied level
 */
Level& Level::operator=(const Level &rhs)
{
  // Operator= should physically copy pixels, not use CVD's reference-counting image copy.
  ImageRef rhsSize = rhs.im.size();
  if (rhsSize.x != 0 && rhsSize.y != 0 )
  {
    im.resize(rhs.im.size());
    copy(rhs.im, im);
  }

  vCorners = rhs.vCorners;
  vMaxCorners = rhs.vMaxCorners;
  vCornerRowLUT = rhs.vCornerRowLUT;
  return *this;
}

// -------------------------------------------------------------
// Some useful globals defined in LevelHelpers.h live here:
Vector<3> gavLevelColors[LEVELS];

// These globals are filled in here. A single static instance of this struct is run before main()
struct LevelHelpersFiller // Code which should be initialised on init goes here; this runs before main()
{
  LevelHelpersFiller()
  {
    for(int i=0; i<LEVELS; i++)
    {
      if(i==0)  gavLevelColors[i] = makeVector( 1.0, 0.0, 0.0);
      else if(i==1)  gavLevelColors[i] = makeVector( 1.0, 1.0, 0.0);
      else if(i==2)  gavLevelColors[i] = makeVector( 0.0, 1.0, 0.0);
      else if(i==3)  gavLevelColors[i] = makeVector( 0.0, 0.0, 0.7);
      else gavLevelColors[i] =  makeVector( 1.0, 1.0, 0.7); // In case I ever run with LEVELS > 4
    }
  }
};
static LevelHelpersFiller foo;
