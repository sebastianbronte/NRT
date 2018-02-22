/* 
 * File:   DescriptorMatcher.h
 * Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
 *
 * Created on Oct 2016
 * 
 * This class implements the Descriptor based matching based on OpenCV
 * 
 */

#ifndef DESCRIPTORMATCHER_OCV_H
#define	DESCRIPTORMATCHER_OCV_H

//opencv includes
#include <opencv/cxcore.h>
#include <opencv/cv.h>
#include <opencv2/features2d.hpp>

class DescriptorMatcherOCV
{
  public:

    enum Feature2DType{
      ORB,
      BRISK,
      MSER,
      KAZE,
      AKAZE,
      SIFT,
      SURF
    };

    enum MatchType{
      BRUTEFORCE,
      BRUTEFORCE_L1,
      BRUTEFORCE_SL2,
      BRUTEFORCE_HAMMING,
      BRUTEFORCE_HAMMING2,
      FLANNBASED
    };

    DescriptorMatcherOCV(bool debug);
    ~DescriptorMatcherOCV();

    void extractFeaturesAndDescriptorsFromImage(cv::Mat &inImgCurr);

    //thread safe getters for kpts and descriptors.
    std::vector<cv::KeyPoint> & getPoints();
    cv::Mat & getDescriptor();

    void match(cv::Mat &d1, 
               std::vector<cv::KeyPoint> &kp1,
               cv::Mat &d2,
               std::vector<cv::KeyPoint> &kp2,
               std::vector<cv::DMatch> &mtch,
               float distLimit);

    unsigned int getDescriptorLength(){return descriptor->descriptorSize();}
    bool isBinaryType();
    bool isFloatType();

  private:

    cv::Ptr<cv::Feature2D> descriptor;
    //TODO: in case it is a MSER descriptor, it needs a CvMSERParam instead of a feature2D matcher ...
    cv::Ptr<cv::DescriptorMatcher> matcher;
    Feature2DType featureType;
    MatchType matchType;

    std::vector<cv::KeyPoint> m_kpts;

    cv::Mat m_desc;

    bool mbDebug;

    void descriptorTypeParser();
    void matchingTypeParser();

    //TODO: add LUT/LUTs to speed up the finding process
    void checkAndAddMatch(std::vector<cv::DMatch> &mtch, const cv::DMatch &math);
};

#endif	/* DESCRIPTORMATCHER_OCV_H */
