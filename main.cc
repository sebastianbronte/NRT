// Copyright 2008 Isis Innovation Limited
// This is the main extry point for NRT
#include <stdlib.h>
#include <iostream>
#include <gvars3/instances.h>
#include "System.h"

using namespace std;
using namespace GVars3;

int main(int argc, char *argv[])
{
  cout << "  Welcome to NRT " << endl;
  cout << "  --------------- " << endl;
  cout << "  Real-time Non-Rigid Tracking" << endl;
  cout << endl;
  cout << "  Parsing settings.cfg ...." << endl;
  GUI.LoadFile("settings.cfg");
  cout << "  done" << endl;

  //Checking input and output data and video streams

  cout << "  * Checking video input" << endl;
  VideoSource *mpVideoSource = NULL;
  string inputopt = GV3::get<string>("i","cam");
  string inputseq;
  try{
      if(!inputopt.compare("cam")){
        mpVideoSource = new VideoSource_DV();
        cout << "    - Live camera video input is selected" << endl;
      }else if(!inputopt.compare("wcam")){
        mpVideoSource = new VideoSource_V4L();
        cout << "    - Live webcam video input is selected" << endl;
      }else if(!inputopt.compare("rec")){
        string opt2 = GV3::get<string>("ifile");
        mpVideoSource = new VideoSource_FF(opt2);
        cout << "    - Recorded input sequence is selected: " << opt2 << endl;
        string opt3 = GV3::get<string>("coeffs","");
        if(!opt3.empty())
            dynamic_cast<VideoSource_FF*>(mpVideoSource)->GetCoefficientsFromFile(opt3);
      }else if(!inputopt.compare("txt")){
        string opt2 = GV3::get<string>("tracks","");
        if(!opt2.empty())
        {
            cout << "    - Plain text file with input tracks is selected" << endl;
            mpVideoSource = new VideoSource_TXT(opt2);
            string opt3 = GV3::get<string>("visfile","");
            if(!opt3.empty())
                dynamic_cast<VideoSource_TXT*>(mpVideoSource)->GetVisibilityFromFile(opt3);
            opt3 = GV3::get<string>("coeffs","");
            if(!opt3.empty())
                dynamic_cast<VideoSource_TXT*>(mpVideoSource)->GetCoefficientsFromFile(opt3);
        }
        else
        {
            cout << " A valid text filename for the tracks must be given. exiting" << endl;
            return 0;
        }
      }else{
        cout << "  A valid filename must be given. exiting" << endl;
        return 0;
      }
  }catch(CVD::Exceptions::All e){
      cout << endl;
      cout << "!! Failed when initializing the system; got exception. " << endl;
      cout << "   Exception was: " << endl;
      cout << e.what << endl;
  }

  GUI.StartParserThread(); // Start parsing of the console input
  atexit(GUI.StopParserThread); 

  try
  {
      System s(mpVideoSource);
      s.Run();
  }
  catch(CVD::Exceptions::All e)
  {
      cout << endl;
      cout << "!! Failed to run system; got exception. " << endl;
      cout << "   Exception was: " << endl;
      cout << e.what << endl;
  }
  return 0;
}
