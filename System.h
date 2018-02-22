// -*- c++ -*-
// Copyright 2008 Isis Innovation Limited
//
// Modified by Sebastian Bronte <sebastian.bronte@depeca.uah.es>
//
// System.h
//
// Defines the System class
//
// This stores the main functional classes of the system, like the
// mapmaker, map, tracker etc, and spawns the working threads.
//
#ifndef __SYSTEM_H
#define __SYSTEM_H
#include "VideoSource.h"
#include "GLWindow2.h"
#include "DBases.h"

#include <gvars3/instances.h>

#include <cvd/image.h>
#include <cvd/rgb.h>
#include <cvd/byte.h>

class ATANCamera;
class Map;
class AbstractTracker;
//class VisualNRTracker;
//class NRTracker;
//class Tracker;
//class ARDriver;
class MapViewer;
class TrackWriter;

class System
{
public:
  System(VideoSource * mpVideoSource);
  ~System();

  void Run();

private:
  VideoSource *mVideoSource; //converting this member to a pointer in order to be able to change the input stream by means of polymorphism
  GLWindow2 mGLWindow;
  CVD::Image<CVD::Rgb<CVD::byte> > mimFrameRGB;
  CVD::Image<CVD::byte> mimFrameBW;

  Map *mpMap;
//  VisualNRTracker_v1 *mpTracker;
//  VisualNRTracker *mpTracker;
//  NRTracker *mpTracker;
  AbstractTracker *mpTracker;
  ATANCamera *mpCamera;
  //ARDriver *mpARDriver;
  MapViewer *mpMapViewer;
  TrackWriter *mpTrackWriter;
  DBases *mpBase;

#ifdef _LINUX
    GVars3::gvar3<int> mgvnSaveFIFO;                // Output to a FIFO (make a video)
    GVars3::gvar3<int> mgvnBitrate;                 // Bitrate to encode at
#endif

  bool mbDone;

  static void GUICommandCallBack(void* ptr, std::string sCommand, std::string sParams);
  void SaveFIFO();                                // save the video out to a FIFO (save to disk)
};



#endif
