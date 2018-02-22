// -*- c++ -*-
// Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
// Date: 5/9/13
//
// DeformablePatch.h
//
// Declares DeformablePatch class
// 
// This is a more ellaborated hierarchical patch matching class,
// used when basic PatchFinder does not work properly. tracking small patches
// it's used by the tracker for building the initial map

#ifndef __DEFORMABLE_PATCH_H
#define __DEFORMABLE_PATCH_H

#include <cvd/image.h>
#include <cvd/byte.h>
#include <cvd/utility.h>
#include <TooN/TooN.h>
using namespace TooN;
#include <vector>

#include "MiniPatch.h"
#include "PatchFinder.h"

class KeyFrame;

struct DeformablePatch
{
  public:
    DeformablePatch():mFinder(){};

    // Finds patch in a new image
    bool FindPatch(CVD::ImageRef &irPos,
                   KeyFrame &kf,
                   int nRange);

    // Second version to look for a patch.
    bool FindPatch2(CVD::ImageRef &ir, KeyFrame &kf1, KeyFrame &kf2, int nRange, bool showStats=false);
    //Spiral search version
    bool FindPatch3(CVD::ImageRef &ir1, CVD::ImageRef &ir2, KeyFrame &kf1, KeyFrame &kf2, int nRange);

   CVD::ImageRef getLastPos(){return mirLastFound;}
   int getLastLevel(){return mnLastLevel;}

  protected:
  
    void preFilterAreaCorners(KeyFrame &kf,
                              std::vector<CVD::ImageRef>& corners,
                              CVD::ImageRef &irPos,
                              int nRange,
                              int minRadius);
    bool basicSearch(CVD::ImageRef &irPos,
                     CVD::ImageRef &irFound,
                     KeyFrame &kf1,
                     KeyFrame &kf2,
                     int nRange,
                     std::vector<CVD::ImageRef>&filteredCorners);

    MiniPatch mFinder;

    CVD::ImageRef mirLastFound;
    int mnLastLevel;
};

#endif
