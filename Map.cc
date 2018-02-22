// Copyright 2008 Isis Innovation Limited
#include "Map.h"
#include "NRTrackerData.h"
#include "NRMapPoint.h"

/**
 * Constructor. Calls reset and sets the map ID
 */
Map::Map()
{
  Reset();
}

Map::~Map()
{
  Reset();
}
/**
 * Reset the map
 */
void Map::Reset()
{
  for(unsigned int i=0; i<vpPoints.size(); i++)
  {
    delete vpPoints[i]->pTData;
    delete vpPoints[i];
  }
  vpPoints.clear();
  bGood = false;
  EmptyTrash();
}

/**
 * Move any points marked as bad to the trash
 */
void Map::MoveBadPointsToTrash()
{
  int nBad = 0;
  for(int i = vpPoints.size()-1; i>=0; i--)
  {
    if(vpPoints[i]->bBad)
    {
      vpPointsTrash.push_back(vpPoints[i]);
      vpPoints.erase(vpPoints.begin() + i);
      nBad++;
    }
  }
}

/**
 * Delete of the points in the trash
 */
void Map::EmptyTrash()
{
  for(unsigned int i=0; i<vpPointsTrash.size(); i++)
  {
    delete vpPointsTrash[i]->pTData;
    delete vpPointsTrash[i];
  }
  vpPointsTrash.clear();
}
