// -*- c++ *--
// Copyright 2008 Isis Innovation Limited
//
// VideoSource.h
// Declares the VideoSource class
// 
// This is a very simple class to provide video input; this can be
// replaced with whatever form of video input that is needed.  It
// should open the video input on construction, and provide two
// function calls after construction: Size() must return the video
// format as an ImageRef, and GetAndFillFrameBWandRGB should wait for
// a new frame and then overwrite the passed-as-reference images with
// GreyScale and Colour versions of the new frame.
//

// Modified by Sebastian Bronte to allow multiple VideoSource on configuration,
// instead of compilation time.

#define __STDC_CONSTANT_MACROS //to make Ffmpeg stuff compile correctly

#include <cvd/image.h>
#include <cvd/byte.h>
#include <cvd/rgb.h>
#include <TooN/TooN.h>

struct VideoSourceData;

class VideoSource
{
  public:
    VideoSource():mptr(NULL){};
    virtual ~VideoSource(){}
    virtual void GetAndFillFrameBWandRGB(CVD::Image<CVD::byte> &imBW, CVD::Image<CVD::Rgb<CVD::byte> > &imRGB)=0;
    virtual void GetCoefficientsCurrentFrame(unsigned int nframe, std::vector<double> &coeffs)=0;
    virtual CVD::ImageRef Size()=0;

  protected:
    void *mptr;
    CVD::ImageRef mirSize;
};

class VideoSource_FF: public VideoSource
{
  public:
    VideoSource_FF(std::string filename);
    virtual ~VideoSource_FF();
    void GetAndFillFrameBWandRGB(CVD::Image<CVD::byte> &imBW, CVD::Image<CVD::Rgb<CVD::byte> > &imRGB);
    void GetCoefficientsFromFile(std::string filename);
    void GetCoefficientsCurrentFrame(unsigned int nframe, std::vector<double> &coeffs);
    CVD::ImageRef Size();
  private:
    unsigned int nFrames;
    unsigned int nCoeffs;
    double *mpcoeffs;
};

class VideoSource_DV: public VideoSource
{
  public:
    VideoSource_DV();
    virtual ~VideoSource_DV();
    void GetAndFillFrameBWandRGB(CVD::Image<CVD::byte> &imBW, CVD::Image<CVD::Rgb<CVD::byte> > &imRGB);
    void GetCoefficientsCurrentFrame(unsigned int nframe, std::vector<double> &coeffs){}
    CVD::ImageRef Size();
};

class VideoSource_V4L: public VideoSource
{
  public:
    VideoSource_V4L();
    virtual ~VideoSource_V4L();
    void GetAndFillFrameBWandRGB(CVD::Image<CVD::byte> &imBW, CVD::Image<CVD::Rgb<CVD::byte> > &imRGB);
    void GetCoefficientsCurrentFrame(unsigned int nframe, std::vector<double> &coeffs){}
    CVD::ImageRef Size();
};

class VideoSource_TXT: public VideoSource
{
  public:
    VideoSource_TXT(std::string filename);
    ~VideoSource_TXT();
    void GetAndFillFrameBWandRGB(CVD::Image<CVD::byte> &imBW, CVD::Image<CVD::Rgb<CVD::byte> > &imRGB);
    CVD::ImageRef Size();
    void GetTracksCurrentFrame(unsigned int nframe, std::vector<TooN::Vector <2> > &tracks);
    void GetVisibilityCurrentFrame(unsigned int nframe, std::vector<bool> &visflags);
    void GetVisibilityFromFile(std::string filename);
    void GetCoefficientsFromFile(std::string filename);
    void GetCoefficientsCurrentFrame(unsigned int nframe, std::vector<double> &coeffs);
    inline unsigned int GetNFrames(){return nFrames;}
    inline unsigned int GetNPoints(){return nPoints;}
    inline bool areCoefficientsLoaded(){return mpcoeffs==NULL;}
  private:
    bool GetTracksFromFile(std::string filename);
    unsigned int nFrames;
    unsigned int nPoints;
    unsigned int nCoeffs;
    bool *mpvis;
    double *mpcoeffs;
};
