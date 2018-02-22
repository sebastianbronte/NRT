// Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>

#include "PatchFinder.h"

// Date: 5/9/13

#include "DeformablePatch.h"
#include "KeyFrame.h"
#include <gvars3/instances.h>
#include <cvd/fast_corner.h>

using namespace CVD;
using namespace std;
using namespace GVars3;

// Find a patch by searching at FAST corners in an input image
// If available, a row-corner LUT is used to speed up search through the
// FAST corners
bool DeformablePatch::FindPatch(ImageRef &irPos,
			  KeyFrame &kf,
			  int nRange)
{
  ImageRef irPosCenter = irPos;

  if (!kf.aLevels[0].im.in_image_with_border(irPos, mFinder.mnHalfPatchSize))
      return false;

  ImageRef disp13(nRange/4, nRange/4);
  ImageRef disp24(-nRange/4, nRange/4);

  //  ______________
  //  |      |      |
  //  |  p2  | p1   |
  //  |------|------|
  //  |  p3  | p4   |
  //  |______|______|

  ImageRef irPos1 = irPosCenter + disp13;
  ImageRef irPos1_2 = irPos1;
  ImageRef irPos2 = irPosCenter + disp24;
  ImageRef irPos2_2 = irPos2;
  ImageRef irPos3 = irPosCenter - disp13;
  ImageRef irPos3_2 = irPos3;
  ImageRef irPos4 = irPosCenter - disp24;
  ImageRef irPos4_2 = irPos4;

  //prefilter corners and decide which matching function to use

  mFinder.SampleFromImage(irPos1,kf.aLevels[0].im);
  bool found1 = mFinder.FindPatch(irPos1, kf.aLevels[0].im, nRange/2,kf.aLevels[0].vCorners);
  mFinder.SampleFromImage(irPos2,kf.aLevels[0].im);
  bool found2 = mFinder.FindPatch(irPos2, kf.aLevels[0].im, nRange/2,kf.aLevels[0].vCorners);
  mFinder.SampleFromImage(irPos3,kf.aLevels[0].im);
  bool found3 = mFinder.FindPatch(irPos3, kf.aLevels[0].im, nRange/2,kf.aLevels[0].vCorners);
  mFinder.SampleFromImage(irPos4,kf.aLevels[0].im);
  bool found4 = mFinder.FindPatch(irPos4, kf.aLevels[0].im, nRange/2,kf.aLevels[0].vCorners);

  cout << "central Pos " << irPos << endl;
  cout << "upper right: " << (found1?"":"not") << " found | irPosFound " << irPos1 << " irPos " << irPos1_2 << endl;
  cout << "upper left: " << (found2?"":"not") << " found | irPosFound " << irPos2 << " irPos " << irPos2_2 << endl;
  cout << "lower left: " << (found3?"":"not") << " found | irPosFound " << irPos3 << " irPos " << irPos3_2 << endl;
  cout << "lower right: " << (found4?"":"not") << " found | irPosFound " << irPos4 << " irPos " << irPos4_2 << endl;

  return found1 || found2 || found3 || found4;
}

bool DeformablePatch::FindPatch2(ImageRef &ir, KeyFrame &kf1, KeyFrame &kf2, int nRange, bool showStats){
    ImageRef irPosCenter = ir;

    if(!kf1.aLevels[0].im.in_image_with_border(ir, nRange/2))
      return false;

    bool found = false;
    mnLastLevel = -1;
    mirLastFound = ir;

    //first step. Lowest level.
    int l = kf1.MAXLEVELS - 1;

    for(l=kf1.MAXLEVELS-1;l>=0;l--)
    {
      int nLevelScale = LevelScale(l);
      irPosCenter = mirLastFound / nLevelScale;
      int nRangel = 8;
      if(l == kf1.MAXLEVELS -1){
        nRangel = (nRange + nLevelScale - 1) / nLevelScale;
      }else{
        nRangel = nRange / (nLevelScale);
        nRangel = (nRangel<15)?nRangel:15;
      }
      
      nRangel = nRangel<4?4:nRangel;

      if(!kf2.aLevels[l].im.in_image_with_border(irPosCenter, nRangel) ||
         !kf1.aLevels[l].im.in_image_with_border(irPosCenter, nRangel))
        continue;

      mFinder.SampleFromImage(irPosCenter, kf1.aLevels[l].im);
      ImageRef irPosCenter2 = irPosCenter;
      //pure correlation matching
      found = mFinder.FindPatchNoCorners(irPosCenter2, kf2.aLevels[l].im, nRangel, showStats);

      if( l == (kf1.MAXLEVELS -1)){
        if(found){
          irPosCenter = irPosCenter2;
          mirLastFound = irPosCenter * nLevelScale;
          mnLastLevel = l;
        }else{
          mirLastFound = ir;
          mnLastLevel = -1;
          break;
        }
      }else{
        if(found){
          irPosCenter = irPosCenter2;
          mirLastFound = irPosCenter * nLevelScale;
          mnLastLevel = l;
        }else{
          mirLastFound = irPosCenter * nLevelScale;
          if(irPosCenter2.x == -1 && irPosCenter2.y == -1)
              mnLastLevel=-1;
        }
      }
    } 
    return mnLastLevel!=-1;
}

bool DeformablePatch::FindPatch3(ImageRef &ir1, ImageRef &ir2, KeyFrame &kf1, KeyFrame &kf2, int nRange)
{
  //improve performance: look for the point from center to the corners of the patch
  ImageRef irPos = ir1, irPos2=ir2;
  if(!kf1.aLevels[0].im.in_image_with_border(irPos,nRange))
  {
      return false;
  }

  vector<ImageRef> vCornersFiltered1, vCornersFiltered2;
  
  preFilterAreaCorners(kf1, vCornersFiltered1, irPos, nRange, 0);
  preFilterAreaCorners(kf2, vCornersFiltered2, irPos, nRange, 0);
  
  bool f=basicSearch(irPos, irPos2, kf1, kf2, nRange, vCornersFiltered2);
  if(f)
  {
    ir1 = irPos;
    ir2 = irPos2;
    return true;
  }
  
  unsigned int v=0;
  vector<unsigned int> vCornerRowLUT; vCornerRowLUT.clear();
  for(int y=irPos.y-nRange; y<irPos.y + nRange; y++)
  {
    while( (v < vCornersFiltered1.size()) && (y > vCornersFiltered1[v].y) )
      v++;
    vCornerRowLUT.push_back(v);
  }
  
  //spiral reading of the matrix from irPos to find correspondences
  int r=1;
  vector<ImageRef>::iterator it;
  do
  {
    //left to right
    int rowPoints = vCornerRowLUT[nRange-r] - vCornerRowLUT[nRange-r-1];
    if( rowPoints > 0)
    {
      for(int p=0;p<rowPoints && !f;p++)
      {
        it = vCornersFiltered1.begin() + vCornerRowLUT[nRange-r] + p;
        //in the interval, look for the point
        if(it->x >= irPos.x-r && it->x <= irPos.x+r)
        {
          irPos2 = *it;
          irPos = irPos2;
          f = basicSearch(*it, irPos2, kf1, kf2, nRange, vCornersFiltered2);
        }
      }
    }
    if(f) break;
    //up to down
    rowPoints = vCornerRowLUT[nRange+r] - vCornerRowLUT[nRange-r];
    if( rowPoints > 0)
    {
      for(int p = 0;p < rowPoints && !f; p++)
      {
        it = vCornersFiltered1.begin() + vCornerRowLUT[nRange-r] + p;
        //in the interval, look for the point
        if(it->y <= irPos.y+r && it->y >= irPos.y-r && (it->x == irPos.x-r || it->x == irPos.x+r))
        {
          irPos2 = *it;
          irPos = irPos2;
          f = basicSearch(*it, irPos2, kf1, kf2, nRange, vCornersFiltered2);
        }
      }
    }
    if(f) break;
    //right to left
    rowPoints = vCornerRowLUT[nRange+r] - vCornerRowLUT[nRange+r-1];
    if( rowPoints > 0)
    {
      for(int p=0;p<rowPoints && !f;p++)
      {
        it = vCornersFiltered1.begin() + vCornerRowLUT[nRange+r] + p;
        //in the interval, look for the point
        if(it->x >= irPos.x-r && it->x <= irPos.x+r) 
        {
          irPos2 = *it;
          irPos = irPos2;
          f = basicSearch(*it, irPos2, kf1, kf2, nRange, vCornersFiltered2);
        }
      }
    }
    if(f) break;
    r++;       
  }while(r<nRange);
  if(f)
  {
    ir1 = irPos;
    ir2 = irPos2;
    return true;
  }  
  return false;
}

bool DeformablePatch::basicSearch(ImageRef &irPos, ImageRef &irFound, KeyFrame &kf1, KeyFrame &kf2, int nRange, vector<ImageRef>&filteredCornersSecond)
{
  if(!kf1.aLevels[0].im.in_image_with_border(irPos,nRange))
  {
      return false;
  }
  mFinder.SampleFromImage(irPos,kf1.aLevels[0].im);
  irFound = irPos;
  return mFinder.FindPatch(irFound, kf2.aLevels[0].im, nRange, filteredCornersSecond);
}

void DeformablePatch::preFilterAreaCorners(KeyFrame &kf, vector<ImageRef>& corners,
                                           ImageRef &irPos, int nRange, int minRadius)
{  
  ImageRef irBest;
  ImageRef irBBoxTL = irPos - ImageRef(nRange, nRange);
  ImageRef irBBoxBR = irPos + ImageRef(nRange, nRange);
  vector<ImageRef>::iterator i;

  for(i = kf.aLevels[0].vCorners.begin(); i!=kf.aLevels[0].vCorners.end(); i++)
    if(i->y >= irBBoxTL.y) break;
  
  for(; i!=kf.aLevels[0].vCorners.end(); i++)
  {
    if(i->x < irBBoxTL.x  || i->x > irBBoxBR.x)
      continue;
    if(i->y > irBBoxBR.y)
      break;
    
    if(minRadius>0)
    {
      //now check against the ones already in corners vector if they are within "gap" distance
      bool f = false;
      for(vector<ImageRef>::iterator j = corners.begin(); j!= corners.end();j++)
        if((j->y-i->y)<minRadius && (j->x-i->x)<minRadius)
	  f=true;
      
      if(!f)
        corners.push_back(*i);
    }else
      corners.push_back(*i);
  }
  
}
