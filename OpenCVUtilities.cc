/*
* Adapter to OpenCV class
* Original idea from https://github.com/BeLioN-github/PTAM/blob/master/VideoSource_Linux_OpenCV.cc
* Arnaud GROSJEAN (VIDE SARL)
* 
* Modifications by Sebastian Bronte for conversions between formats.
* <sebastian.bronte@depeca.uah.es>
*
* INSTALLATION :
* - Copy the VideoSource_Linux_OpenCV.cc file in your PTAM directory
* - In the Makefile:
*	- set the linkflags to
	LINKFLAGS = -L MY_CUSTOM_LINK_PATH -lblas -llapack -lGVars3 -lcvd -lcv -lcxcore -lhighgui
*	- set the videosource to 
	VIDEOSOURCE = VideoSource_Linux_OpenCV.o
* - Compile the project
* - Enjoy !
*/
#include "OpenCVUtilities.h"
#include <cvd/colourspace_convert.h>
#include <cvd/colourspaces.h>
#include <gvars3/instances.h>
#include <opencv/highgui.h>
#include <opencv/cxcore.h>
#include <opencv/cv.h>

OpenCVUtilities::OpenCVUtilities()
{
}

OpenCVUtilities::~OpenCVUtilities()
{
}

/**
 * function to convert a BW image from OpenCV to libcvd format. Not tested.
 * @param frame
 * @param imBW
 */
void OpenCVUtilities::conversionOpenCV2CVD_NB(const cv::Mat &frame, CVD::BasicImage<CVD::byte> &imBW)
{
  int nchannels = frame.dims;
  for (int i = 0; i < frame.rows; i++)
  {
    for (int j = 0; j < frame.cols; j++)
    {
      if(nchannels == 2)
      {
        imBW[i][j] = *frame.ptr<unsigned char>(i,j);
      }
      else
      {
        imBW[i][j] = ((*(frame.ptr<unsigned char>(i,j,0)))+
                      (*(frame.ptr<unsigned char>(i,j,1)))+
                      (*(frame.ptr<unsigned char>(i,j,2))))/3;
      }
    }
  }

}

/**
 * function to convert a RGB image from OpenCV to libcvd format. Not tested
 * @param frame
 * @param imRGB
 */
void OpenCVUtilities::conversionOpenCV2CVD_RGB(const cv::Mat &frame, CVD::BasicImage<CVD::Rgb<CVD::byte> > &imRGB)
{
  cv::Mat_<cv::Vec3b>& frame_p = (cv::Mat_<cv::Vec3b>&)frame;
  for (int i = 0; i < frame.rows; i++)
  {
    for (int j = 0; j < frame.cols; j++)
    {
      imRGB[i][j].red = frame_p(i,j)[2];
      imRGB[i][j].green = frame_p(i,j)[1];
      imRGB[i][j].blue = frame_p(i,j)[0];
    }
  }
}

/**
 * function to convert a BW image from libcvd to OpenCV format. Tested
 * @param imBW
 * @param frame
 */
void OpenCVUtilities::conversionCVD2OpenCV_NB(const CVD::BasicImage<CVD::byte> &imBW, cv::Mat &frame)
{
  CVD::ImageRef size = imBW.size();
  for (int i = 0; i < size.y; i++)
  {
    for(int j = 0; j < size.x; j++)
    {
      frame.at<unsigned char>(i,j) = imBW[i][j];
    }
  }
}

/**
 * function to convert from RGB image in libcvd to OpenCV. Not tested
 * @param imRGB
 * @param frame
 */
void OpenCVUtilities::conversionCVD2OpenCV_RGB(const CVD::BasicImage<CVD::Rgb<CVD::byte> > &imRGB, cv::Mat &frame)
{
  CVD::ImageRef size = imRGB.size();
  for (int i=0;i<size.y;i++)
  {
    for(int j=0;j<size.x;j++)
    {
        frame.at<unsigned char>(i,j,0) = imRGB[i][j].red;
        frame.at<unsigned char>(i,j,1) = imRGB[i][j].green;
        frame.at<unsigned char>(i,j,2) = imRGB[i][j].blue;
    }
  }
}
