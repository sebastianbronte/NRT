/* 
 * File:   OpenCVUtilities.h
 * Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
 *
 * Created on 7 Feb 2014
 * 
 * This class is an adapter to go from OpenCV structures to libCVD structures 
 * and viceversa
 * 
 */

#ifndef OPENCVUTILITIES_H
#define	OPENCVUTILITIES_H

//libcvd includes
#include <cvd/image.h>
#include <cvd/byte.h>
#include <cvd/utility.h>
//toon includes
#include <TooN/TooN.h>
//opencv includes
#include <opencv/cxcore.h>
#include <opencv/cv.h>

class OpenCVUtilities
{
  private:
    OpenCVUtilities();
    ~OpenCVUtilities();
  public:
    //OpenCV to libCVD conversions
    static void conversionOpenCV2CVD_NB(const cv::Mat &frame,
                                        CVD::BasicImage<CVD::byte> &imBW);

    static void conversionOpenCV2CVD_RGB(const cv::Mat &frame,
                                         CVD::BasicImage<CVD::Rgb<CVD::byte> > &imRGB);

    //libCVD to OpenCV conversions
    static void conversionCVD2OpenCV_NB(const CVD::BasicImage<CVD::byte> &imBW,
                                        cv::Mat &frame);

    static void conversionCVD2OpenCV_RGB(const CVD::BasicImage<CVD::Rgb<CVD::byte> > &imRGB,
                                         cv::Mat &frame);

};

#endif	/* OPENCVUTILITIES_H */
