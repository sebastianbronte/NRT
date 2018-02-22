// Copyright 2008 Isis Innovation Limited
#include "System.h"
#include "OpenGL.h"
#include <gvars3/instances.h>
#include <stdlib.h>
#include "ATANCamera.h"
//#include "MapMaker.h"
//#include "NRMapMaker.h"
//#include "NRMapMaker.h"
#include "VisualNRTrackerV1.h"
#include "VisualNRTracker.h"
#include "NRTracker.h"
#include "Tracker.h"
//#include "ARDriver.h"
#include "MapViewer.h"
#include "TrackWriter.h"
//#include "NRTracker.h"
#include "Interpolation.h"
#ifdef _LINUX
#include <fcntl.h>
#endif
#include<sys/time.h>

using namespace CVD;
using namespace std;
using namespace GVars3;


System::System(VideoSource * mpVideoSource)
  : mVideoSource(mpVideoSource),mGLWindow(mVideoSource ? mVideoSource->Size() : ImageRef(640,480), "PTAM")
{
  //check if mpVideoSource is NULL
  ImageRef s;
  if(mVideoSource)
    s=mVideoSource->Size();
  else
    s = ImageRef(640,480);

  mimFrameBW.resize(s);
  mimFrameRGB.resize(s);

  GUI.RegisterCommand("exit", GUICommandCallBack, this);
  GUI.RegisterCommand("quit", GUICommandCallBack, this);
  
  // First, check if the camera is calibrated.
  // If not, we need to run the calibration widget.
  Vector<NUMTRACKERCAMPARAMETERS> vTest;

#ifdef _LINUX
  GV2.Register(mgvnSaveFIFO, "SaveFIFO", 0, SILENT);
  GV2.Register(mgvnBitrate, "Bitrate", 15000, SILENT);
#endif

  vTest = GV3::get<Vector<NUMTRACKERCAMPARAMETERS> >("Camera.Parameters", ATANCamera::mvDefaultParams, HIDDEN);
  mpCamera = new ATANCamera("Camera");
  
  if(vTest == ATANCamera::mvDefaultParams)
  {
    cout << endl;
    cout << "! Camera.Parameters is not set, need to run the CameraCalibrator tool" << endl;
    cout << "  and/or put the Camera.Parameters= line into the appropriate .cfg file." << endl;
    exit(1);
  }
  mpBase = new DBases();
  string basedisk = GV3::get<string>("bases","");

  if(!basedisk.empty())
  {
    if (!mpBase->loadFromFile(basedisk.c_str()))
    {
      Exceptions::All e;
      cout << "error while trying to load the bases file" << endl;
      delete mpBase;
      exit(1);
    }
    cout << "Bases successfully loaded from file" << endl;

    string rigid = GV3::get<string>("rigid","");
    if(!rigid.empty())
    {
      if(!mpBase->loadRigidFromFile(rigid))
      {
          Exceptions::All e;
          cout << "error while loading the rigid initial shape file" << endl;
          delete mpBase;
          exit(1);
      }
      cout << "Rigid shape file loaded successfully" << endl;
    }

  }
  mpMap = new Map;

  string trackingType = GV3::get<string>("Tracker.Type","");
  if (trackingType == "Points")
    mpTracker = new NRTracker(s, *mpCamera, *mpMap, *mpBase);
  else if (trackingType == "Descriptor")
    mpTracker = new VisualNRTracker(s, *mpCamera, *mpMap, *mpBase);
  else if (trackingType == "PTAM-based")
    mpTracker = new VisualNRTrackerV1(s, *mpCamera, *mpMap, *mpBase);
  else
    mpTracker = new Tracker(s, *mpCamera, *mpMap);
  //mpARDriver = new ARDriver(*mpCamera, mVideoSource->Size(), mGLWindow);
  mpMapViewer = new MapViewer(*mpMap, mGLWindow);
  mpTrackWriter = new TrackWriter(*mpTracker);

  GUI.ParseLine("GLWindow.AddMenu Menu Menu");
  GUI.ParseLine("Menu.ShowMenu Root");
  GUI.ParseLine("Menu.AddMenuButton Root Reset Reset Root");
  //GUI.ParseLine("Menu.AddMenuButton Root Spacebar PokeTracker Root");
#ifdef _LINUX
  GUI.ParseLine("MapsMenu.AddMenuToggle Serial \"Save Video\" SaveFIFO Serial");
  GUI.ParseLine("MapsMenu.AddMenuSlider Serial Bitrate Bitrate 100 20000 Serial");
#endif
  //GUI.ParseLine("DrawAR=0");
  GUI.ParseLine("DrawMap=0");
  GUI.ParseLine("Menu.AddMenuToggle Root \"View Map\" DrawMap Root");
  //GUI.ParseLine("Menu.AddMenuToggle Root \"Draw AR\" DrawAR Root");
  mbDone = false;
}

System::~System(){
  if(mpTracker) { delete mpTracker; mpTracker=NULL; }
  //delete mpARDriver;
  if(mpMapViewer) { delete mpMapViewer; mpMapViewer=NULL; }
  if(mpTrackWriter) { delete mpTrackWriter; mpTrackWriter=NULL; }
  if(mpMap) { delete mpMap; mpMap=NULL; }
  if(mpCamera) { delete mpCamera; mpCamera=NULL; }
  if(mpBase) { delete mpBase; mpBase=NULL; }
  if(mVideoSource)
    delete mVideoSource;
}

void System::Run()
{

#ifndef WIN32

  struct sched_param p;
  p.sched_priority = sched_get_priority_max(SCHED_FIFO); // Obtener prioridad maxima
  sched_setscheduler(0, SCHED_FIFO, &p);

#endif

  vector<Vector <2> > tracks(0);
  vector<bool> visibility(0);
  vector<double> coeffs(0);
  char numberstr[7];
  struct timeval t1,t2,tf1,tf2;
  while(!mbDone)
  {
    // We use two versions of each video frame:
    // One black and white (for processing by the tracker etc)
    // and one RGB, for drawing.

    // Grab new video frame...

    gettimeofday(&t1,NULL);

    try
    {

      if(!mpTracker->AreTracksLoadedFromDisk())
        mVideoSource->GetAndFillFrameBWandRGB(mimFrameBW, mimFrameRGB);
      else
      {
        dynamic_cast<VideoSource_TXT*>(mVideoSource)->GetTracksCurrentFrame(mpTracker->GetFrameNumber(), tracks);
        dynamic_cast<VideoSource_TXT*>(mVideoSource)->GetVisibilityCurrentFrame(mpTracker->GetFrameNumber(), visibility);
        dynamic_cast<VideoSource_TXT*>(mVideoSource)->GetCoefficientsCurrentFrame(mpTracker->GetFrameNumber(),coeffs);
      }

      if(mpTracker->AreCoeffsLoadedFromDisk())
      {
         mVideoSource->GetCoefficientsCurrentFrame(mpTracker->GetFrameNumber(), coeffs);
      }
    }
    catch(CVD::Exceptions::All e)
    {
      cout << "End of video sequence" << endl;
      break;
    }

    static bool bFirstFrame = true;
    if(bFirstFrame)
    {
      //mpARDriver->Init();
      bFirstFrame = false;
    }

    mGLWindow.SetupViewport();
    mGLWindow.SetupVideoOrtho();
    mGLWindow.SetupVideoRasterPosAndZoom();

    //if(!mpMap->IsGood())
    //	mpARDriver->Reset();

    static gvar3<int> gvnDrawMap("DrawMap", 0, HIDDEN|SILENT);
    //static gvar3<int> gvnDrawAR("DrawAR", 0, HIDDEN|SILENT);

    bool bDrawMap = mpMap->IsGood() && *gvnDrawMap;
    //bool bDrawAR = mpMap->IsGood() && *gvnDrawAR;
    gettimeofday(&tf1,NULL);
    try
    {
      if(!mpTracker->AreTracksLoadedFromDisk())
      {
        if(!mpTracker->AreCoeffsLoadedFromDisk())
        {
          mpTracker->TrackFrame(mimFrameBW, true/*!bDrawAR && !bDrawMap*/);
        }else{
          mpTracker->TrackFrame(mimFrameBW, coeffs, true/*!bDrawAR && !bDrawMap*/);
        }
      }else{
        mpTracker->TrackFrame(tracks, visibility, coeffs);
      }
    }catch(CVD::Exceptions::All e){
      cout << e.what << endl;
      break;
    }
    gettimeofday(&tf2,NULL);

    if(bDrawMap)
    {
      mpMapViewer->DrawMap(mpTracker->GetCurrentPose());
      if (mpTracker->getInterpolation())
        mpMapViewer->DrawBasesShape(mpBase->mse3RelativePose,
                                    mpTracker->getInterpolation()->m3Shape);
    }
    /*else if(bDrawAR)
      mpARDriver->Render(mimFrameRGB, mpTracker->GetCurrentPose());*/

    //      mGLWindow.GetMousePoseUpdate();
    string sCaption;
    if(bDrawMap)
      sCaption = mpMapViewer->GetMessageForUser();
    else
      sCaption = mpTracker->GetMessageForUser();

    sCaption.append("\n Frame: ");

    sprintf(numberstr,"%05d ", mpTracker->GetFrameNumber());
    sCaption.append(numberstr);

    mGLWindow.DrawCaption(sCaption);
    mGLWindow.DrawMenus();
    mGLWindow.swap_buffers();
    mGLWindow.HandlePendingEvents();

    gettimeofday(&t2,NULL);
    cout << "total frame processing time: " << (t2.tv_sec-t1.tv_sec)*1000 + (t2.tv_usec-t1.tv_usec)*0.001 << " ms" << endl;
    cout << "total tracking processing time: " << (tf2.tv_sec-tf1.tv_sec)*1000 + (tf2.tv_usec-tf1.tv_usec)*0.001 << " ms" << endl;

#ifdef _LINUX
    if( *mgvnSaveFIFO )
    {
      SaveFIFO();
    }
#endif

    }
}

void System::GUICommandCallBack(void *ptr, string sCommand, string sParams)
{
  if(sCommand=="quit" || sCommand == "exit")
    static_cast<System*>(ptr)->mbDone = true;
}

/**
 * Save the current frame to a FIFO.
 * This function is called on each frame to create a video.
 * The GVar SaveFIFO starts and stops the saving, and the GVar
 * Bitrate sets the quality.
 * Bitrate can only be set before the first call of SaveFIFO.
 * @param void
 * @return void
 */
void System::SaveFIFO()
{
#ifdef _LINUX
  //Some static variables
  static CVD::byte* pcImage = NULL;
  static int fd = 0;
  static bool bFIFOInitDone = false;
  static ImageRef irWindowSize;

  if( !bFIFOInitDone )
  {
    irWindowSize = mGLWindow.size();

    ostringstream os;
    os << /*"/bin/bash\n" <<*/
        "file=\"`date '+%Y-%m-%d_%H-%M-%S'`.avi\"; " <<
        "if [ ! -e FIFO ]; then mkfifo FIFO; echo Made FIFO...; fi; " <<
        "echo Mencoding to $file....; " <<
        "cat FIFO |nice mencoder -flip -demuxer rawvideo -rawvideo fps=30:w=" <<
        irWindowSize.x << ":h=" << irWindowSize.y <<
        //":format=rgb24 -o $file -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=" << *mgvnBitrate <<
        ":format=rgb24 -o $file -ovc lavc -lavcopts vcodec=ffv1"
        //":keyint=45 -ofps 30 -ffourcc DIVX - &";
        " -ofps 30 - &";

    int i = system( os.str().c_str() );
    if( i != 0 ) {
      cerr << "ERROR: could not set up the FIFO!" << endl;
      return;
    }

    posix_memalign((void**)(&pcImage), 16, irWindowSize.x*irWindowSize.y*3);
    string s = "FIFO";
    fd = open(s.c_str(), O_RDWR | O_ASYNC);

    bFIFOInitDone = true;
  }

  if( irWindowSize != mGLWindow.size() )
  {
    cerr << "ERROR: Aborting FIFO as window size has changed!!" << endl;
    *mgvnSaveFIFO = 0;
    return;
  }

  glReadBuffer(GL_BACK);
  glReadPixels(0,0,irWindowSize.x,irWindowSize.y,GL_RGB, GL_UNSIGNED_BYTE, pcImage);
  write(fd, (char*) pcImage, irWindowSize.x*irWindowSize.y*3);
#else
  cout << "Video Saving using FIFOs is only available under Linux" << endl;
#endif
}
