/*
 * Author: Sebastian Bronte
 * <sebastian.bronte@depeca.uah.es>
 * Date: Oct 2016
 *
 *  Class that implements the descriptor matching functionality based on OpenCV
 *  It implements the interface to use a variety of descriptors from OpenCV,
 * among others, KAZE, AKAZE, ORB, BRISK, etc.
 */

#include "DescriptorMatcherOCV.h"
#include <cvd/colourspace_convert.h>
#include <cvd/colourspaces.h>
#include <gvars3/instances.h>
#include <iostream>
#include <opencv/cxcore.h>
#include <opencv/cv.h>
#include <opencv2/features2d.hpp>
#include "opencv2/xfeatures2d.hpp"
#include <opencv2/highgui/highgui.hpp>

#include <omp.h>
#include <cvd/thread.h>

#include <string>

using namespace CVD;
using namespace std;
using namespace GVars3;
using namespace cv;

DescriptorMatcherOCV::DescriptorMatcherOCV(bool debug):m_kpts(vector<KeyPoint>()),mbDebug(debug)
{ 
  descriptorTypeParser();
  matchingTypeParser(); 
}

DescriptorMatcherOCV::~DescriptorMatcherOCV()
{
  // Nothing to cleanup in this class
}

void DescriptorMatcherOCV::extractFeaturesAndDescriptorsFromImage(cv::Mat& inImgCurr)
{
  Mat img_32;

  // Convert the images to float
  inImgCurr.convertTo(img_32,CV_32F,1.0/255.0,0);

  //cleanup previous stuff
  std::vector<cv::KeyPoint>().swap(m_kpts);
  m_kpts.clear();
  m_desc.release();

  bool input8bitImg = (featureType==BRISK || featureType==ORB || featureType == SIFT || featureType == SURF);

  if (featureType != MSER){
    descriptor->detectAndCompute(input8bitImg?inImgCurr:img_32, noArray(), m_kpts, m_desc);
  }else{
    //TODO: do stuff to implement MSER?
    cerr << "detection not implemented yet" << endl;
  }

  img_32.release();
}

vector<KeyPoint> & DescriptorMatcherOCV::getPoints()
{
  return m_kpts;
}

Mat & DescriptorMatcherOCV::getDescriptor()
{
  return m_desc;
}

void DescriptorMatcherOCV::match(Mat &d1, vector<KeyPoint> &kp1, 
                                 Mat &d2, vector<KeyPoint> &kp2, 
                                 vector<DMatch> &mtch, float distLimit)
{

  // Matching Descriptors
  vector<vector<DMatch> > dmatches;
  dmatches.clear();
  mtch.clear();

  if (mbDebug){
    cout << "kp1 size " << kp1.size() << "kp2 size " << kp2.size() << endl;
    cout << "d1 size (" << d1.rows << ", " << d1.cols << ") d2 size (" << d2.rows << ", " << d2.cols << ")" << endl;
  }
  if(kp1.size()==0 || kp2.size() == 0 || d1.rows == 0 || d2.rows == 0)
    return;

  matcher->knnMatch(d2,d1,dmatches,2);

  float dist1 = 0.0, dist2 = 0.0;
  // The size of the matches vector is the greatest of the 2 point sets.
  for( size_t i = 0; i < dmatches.size(); i++ )
  {
    DMatch match0 = dmatches[i][0];
    DMatch match1 = dmatches[i][1];
    dist1 = match0.distance;
    dist2 = match1.distance;

    // two criteria: descriptor distance and point distance throught all the candidates of a match.
    Point2f pt01 = kp1[match0.trainIdx].pt;
    Point2f pt02 = kp2[match0.queryIdx].pt;
    Point2f pt11 = kp1[match1.trainIdx].pt;
    Point2f pt12 = kp2[match1.queryIdx].pt;

    Point2f diff0 = pt02 - pt01;
    Point2f diff1 = pt12 - pt11;

    if (mbDebug)
      cout << "match 0 " << match0.trainIdx << ", " << match0.queryIdx
           << " d: " << dist1
           << " pt1 " << pt01 << " pt2 " << pt02 << " dist " << diff0
           << " match 1 " << match1.trainIdx << ", " << match1.queryIdx
           << " d: " << dist2 
           << " pt1 " << pt11 << " pt2 " << pt12 << " dist " << diff1 << endl;

    if (pt01.x > 1 && pt02.x > 1 && pt01.y > 1 && pt02.y > 1 && abs(diff0.x) < distLimit && abs(diff0.y) < distLimit){
      if (mbDebug)
        cout << "checking match 0 ";
      checkAndAddMatch(mtch, match0);
    }
    else if (pt11.x > 1 && pt12.x > 1 && pt11.y > 1 && pt12.y > 1 && abs(diff1.x) < distLimit && abs(diff1.y) < distLimit){
      if (mbDebug)
        cout << "checking match 1 ";
      checkAndAddMatch(mtch, match1);
    }
    else
      if(mbDebug)
        cout << "match discarded" << endl;
  }
}

void DescriptorMatcherOCV::descriptorTypeParser(){
  //default AKAZE
  string featureTypeStr = GV3::get<string>("FeatureType","AKAZE");
  if (!featureTypeStr.compare("KAZE")) {
    // Something
    featureType = KAZE;
    bool extended=false;
    bool upright=false;
    float threshold = 0.001f;
    int nOctaves = 2;
    int nOctaveLayers = 2;
    int diffusivity = KAZE::DIFF_PM_G2;
    descriptor = cv::KAZE::create(extended, upright, threshold, nOctaves,
                                  nOctaveLayers, diffusivity);
  } else if (!featureTypeStr.compare("BRISK")) {
    // Something
    featureType = BRISK;
    int thresh=30;
    int octaves=1;
    float patternScale=1.0f;
    descriptor = cv::BRISK::create(thresh, octaves, patternScale);
    
  } else if (!featureTypeStr.compare("ORB")) {
    // Something
    featureType = ORB;
    int nfeatures=GV2.GetInt("MaxInitialTrails", 1000, SILENT);
    float scaleFactor=1.2f;
    int nlevels=8;
    int edgeThreshold=31;
    int firstLevel=0;
    int WTA_K=2;
    int scoreType=ORB::HARRIS_SCORE;
    int patchSize=31;
    int fastThreshold=20;
    descriptor = cv::ORB::create(nfeatures, scaleFactor, nlevels, edgeThreshold,
        firstLevel, WTA_K, scoreType, patchSize, fastThreshold);
  } else if (!featureTypeStr.compare("MSER")) {
    cerr << "not yet implemented!!!" << endl;
    /*int delta = 5;
    int min_area = 60;
    int max_area = 14400;
    float max_variation=0.25;
    float min_diversity=0.2;
    int max_evolution=200;
    double area_threshold=1.01;
    double min_margin=0.003;
    int edge_blur_size=5;
  
    descriptor = cv::MSER::create(delta, min_area, max_area, max_variation,
                                  min_diversity, max_evolution, area_threshold,
                                  min_margin, edge_blur_size);*/
  }//SURF and SIFT are not available in this implementation of opencv directly. They must be installed as 3rd party software... 
  else if (!featureTypeStr.compare("SURF")) {
    // Something
    featureType = SURF;
    double hessianThreshold=10;
    int nOctaves = 1;
    int nOctaveLayers = 1;
    bool extended = false;
    bool upright = false;
    descriptor = cv::xfeatures2d::SURF::create(hessianThreshold, nOctaves, nOctaveLayers, extended, upright);
  } else if (!featureTypeStr.compare("SIFT")) {
    //Something
    featureType = SIFT;
    descriptor = cv::xfeatures2d::SIFT::create();
  } else {      
    featureType = AKAZE;
    int descriptor_type=cv::AKAZE::DESCRIPTOR_MLDB;
    int descriptor_size = 0;
    int descriptor_channels = 1;
    float threshold = 0.001f;
    int octaves = 1;
    int sublevels = 1;
    int diffusivity = KAZE::DIFF_PM_G2;
    descriptor = cv::AKAZE::create(descriptor_type, descriptor_size, descriptor_channels,
                                   threshold, octaves, sublevels, diffusivity);
  }
  if (mbDebug)
    cout << "feature type" << featureType << endl;
}

void DescriptorMatcherOCV::matchingTypeParser(){
  //default BruteForce-Hamming to be coherent with the default AKAZE.
  string matchingType = GV3::get<string>("MatchingType","BruteForce-Hamming");
  if (!matchingType.compare("BruteForce")) {
    matchType = BRUTEFORCE;
  } else if (!matchingType.compare("BruteForce-L1")){
    matchType = BRUTEFORCE_L1;
  } else if (!matchingType.compare("BruteForce_SL2")){
    matchType = BRUTEFORCE_SL2;
  } else if (!matchingType.compare("BruteForce-Hamming(2)")){
    matchType = BRUTEFORCE_HAMMING2;
  } else if (!matchingType.compare("FlannBased")) {
    matchType = FLANNBASED;
  }
  //default: hamming
  else{
    matchType = BRUTEFORCE_HAMMING;
  }

  matcher = DescriptorMatcher::create(matchingType);
  if (mbDebug)
    cout << "matching type " << matchType << endl;

  // checking incompatibilities
  if (isFloatType() && (matchType == BRUTEFORCE_HAMMING || matchType == BRUTEFORCE_HAMMING2)){
    cerr << "incompatible matcher with feature type" << endl;
  }

}

bool DescriptorMatcherOCV::isBinaryType(){
  return featureType == AKAZE || featureType == BRISK || featureType == ORB;
}

bool DescriptorMatcherOCV::isFloatType(){
  return featureType == KAZE || featureType == SIFT || featureType == SURF;
}

void DescriptorMatcherOCV::checkAndAddMatch(vector<DMatch> &mtch, const DMatch &match)
{
  //TODO: add LUT/LUTs to speed up the finding process
  bool alreadyFound = false;
  unsigned int i = 0;
  for (; i < mtch.size(); i++){
    if (mtch[i].trainIdx == match.trainIdx){
      alreadyFound = true;
      break;
    }
  }

  const DMatch &matchToAdd = match;
  // If already found, check which of the matches is the best one
  if (alreadyFound){
    DMatch &v_match = mtch[i];
    if (mbDebug)
      cout << "Match repeated. (" << i << ": " << v_match.trainIdx << "," << v_match.queryIdx << ") ";
    if (v_match.distance < match.distance){
      const_cast<DMatch&>(matchToAdd) = v_match;
      if (mbDebug)
        cout << "Preserving the previous one" << endl;
    }else{
      if (mbDebug)
        cout << "Substituing for the new one" << endl;
    }
    mtch[i] = matchToAdd;
  }else {
    mtch.push_back(matchToAdd);
    if (mbDebug)
      cout << "Match inserted" << endl;
  }
}

