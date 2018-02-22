// Copyright 2008 Isis Innovation Limited
#include "VideoSource.h"
#include <cvd/Linux/v4l2buffer.h>
#include <cvd/colourspace_convert.h>
#include <cvd/colourspaces.h>
#include <gvars3/instances.h>

using namespace CVD;
using namespace std;
using namespace GVars3;

namespace CVD
{
template<>
struct V4L2_Traits<CVD::yuv422>{
  static const unsigned int pix_code=V4L2_PIX_FMT_YUYV;
};
}

VideoSource_V4L::VideoSource_V4L()
{
  cout << "  VideoSource_Linux: Opening video source..." << endl;
  string QuickCamFile = GV3::get<string>("VideoSource.V4LDevice", "/dev/video0");
  V4L2BufferT<yuv422> *pvb = new V4L2BufferT<yuv422>(QuickCamFile.c_str(), true, V4L2BBMselect);
  mirSize = pvb->size();
  mptr = pvb;
  cout << "  ... got video source." << endl;
}

VideoSource_V4L::~VideoSource_V4L()
{
    delete (V4L2BufferT<yuv422>*) mptr;
}

ImageRef VideoSource_V4L::Size()
{
  return mirSize;
}

void VideoSource_V4L::GetAndFillFrameBWandRGB(Image<byte> &imBW, Image<Rgb<byte> > &imRGB)
{
  V4L2BufferT<yuv422>* pvb = (V4L2BufferT<yuv422>*) mptr;
  V4L2FrameT<yuv422>* pVidFrame = static_cast<V4L2FrameT<yuv422>*>(pvb->get_frame());
  convert_image(*pVidFrame, imBW);
  convert_image(*pVidFrame, imRGB);
  pvb->put_frame(pVidFrame);
}
