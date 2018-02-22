// -*- c++ -*-
// Copyright 2008 Isis Innovation Limited
//
// This header declares the Map class.
// This is pretty light-weight: All it contains is
// a vector of MapPoints and a vector of KeyFrames.
//
// N.b. since I don't do proper thread safety,
// everything is stored as lists of pointers,
// and map points are not erased if they are bad:
// they are moved to the trash list. That way
// old pointers which other threads are using are not 
// invalidated!

#ifndef __MAP_H
#define __MAP_H
#include <vector>
#include <TooN/se3.h>
#include <cvd/image.h>

struct NRMapPoint;
struct KeyFrame;

struct Map
{
  Map();
  ~Map();
  inline bool IsGood() {return bGood;}
  void Reset();

  void MoveBadPointsToTrash();
  void EmptyTrash();

  std::vector<NRMapPoint*> vpPoints;
  std::vector<NRMapPoint*> vpPointsTrash;
  std::vector<KeyFrame*> vpKeyFrames;

  bool bGood;
};

#endif

