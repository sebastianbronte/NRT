// -*- c++ -*-
// Copyright 2008 Isis Innovation Limited
//
// MiniPatch.h
//
// Declares MiniPatch class
// 
// This is a simple pixel-patch class, used for tracking small patches
// it's used by the tracker for building the initial map

#ifndef __MINI_PATCH_H
#define __MINI_PATCH_H

#include <cvd/image.h>
#include <cvd/byte.h>
#include <cvd/utility.h>
#include <TooN/TooN.h>
using namespace TooN;
#include <vector>

struct MiniPatch
{
  // Copy pixels out of source image
  void SampleFromImage(CVD::ImageRef irPos, CVD::BasicImage<CVD::byte> &im);
  // Find patch in a new image
  bool FindPatch(CVD::ImageRef &irPos,
                 CVD::BasicImage<CVD::byte> &im, 
                 int nRange,
                 std::vector<CVD::ImageRef> &vCorners,
                 std::vector<int> *pvRowLUT = NULL);

  // Find patch in a new image without taking into account corners
  bool FindPatchNoCorners(CVD::ImageRef &irPos,
                          CVD::BasicImage<CVD::byte> &im,
                          int nRange,
                          bool showStats=false);

  inline int SSDAtPoint(CVD::BasicImage<CVD::byte> &im, const CVD::ImageRef &ir); // Score function
  static int mnHalfPatchSize;     // How big is the patch?
  static int mnRange;             // How far to search? 
  static int mnMaxSSD;            // Max SSD for matches?
  CVD::Image<CVD::byte> mimOrigPatch;  // Original pixels
};

#endif
