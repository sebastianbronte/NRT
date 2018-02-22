// Copyright 2008 Isis Innovation Limited
#include "MiniPatch.h"
#include <gvars3/instances.h>

using namespace CVD;
using namespace std;
using namespace GVars3;

// Scoring function
inline int MiniPatch::SSDAtPoint(CVD::BasicImage<CVD::byte> &im, const CVD::ImageRef &ir)
{
  if(!im.in_image_with_border(ir, mnHalfPatchSize))
    return mnMaxSSD + 1;
  ImageRef irImgBase = ir - ImageRef(mnHalfPatchSize, mnHalfPatchSize);
  int nRows = mimOrigPatch.size().y;
  int nCols = mimOrigPatch.size().x;
  byte *imagepointer;
  byte *templatepointer;
  int nDiff;
  int nSumSqDiff = 0;
  for(int nRow = 0; nRow < nRows; nRow++)
  {
    imagepointer = &im[irImgBase + ImageRef(0,nRow)];
    templatepointer = &mimOrigPatch[ImageRef(0,nRow)];
    for(int nCol = 0; nCol < nCols; nCol++)
    {
      nDiff = imagepointer[nCol] - templatepointer[nCol];
      nSumSqDiff += nDiff * nDiff;
    }
  }
  return nSumSqDiff;
}

// Find a patch by searching at FAST corners in an input image
// If available, a row-corner LUT is used to speed up search through the
// FAST corners
bool MiniPatch::FindPatch(CVD::ImageRef &irPos, 
                          CVD::BasicImage<CVD::byte> &im, 
                          int nRange, 
                          vector<ImageRef> &vCorners,
                          std::vector<int> *pvRowLUT)
{
  if(!im.in_image_with_border(irPos,nRange))
  {
      return false;
  }
  ImageRef irBest;
  int nBestSSD = mnMaxSSD + 1;
  ImageRef irBBoxTL = irPos - ImageRef(nRange, nRange);
  ImageRef irBBoxBR = irPos + ImageRef(nRange, nRange);
  vector<ImageRef>::iterator i;
  if(!pvRowLUT)
  {
    for(i = vCorners.begin(); i!=vCorners.end(); i++)
      if(i->y >= irBBoxTL.y) break;
  }
  else
  {
    int nTopRow = irBBoxTL.y;
    if(nTopRow < 0)
      nTopRow = 0;
    if(nTopRow >= (int) pvRowLUT->size())
      nTopRow = (int) pvRowLUT->size() - 1;
    i = vCorners.begin() + (*pvRowLUT)[nTopRow];
  }

  for(; i!=vCorners.end(); i++)
  {
    if(i->x < irBBoxTL.x  || i->x > irBBoxBR.x)
      continue;
    if(i->y > irBBoxBR.y)
      break;
    int nSSD = SSDAtPoint(im, *i);

    if(nSSD < nBestSSD)
    {
      irBest = *i;
      nBestSSD = nSSD;
    }
  }

  if(nBestSSD < mnMaxSSD)
  {
    irPos = irBest;
    return true;
  }
  else
    return false;
}

bool MiniPatch::FindPatchNoCorners(CVD::ImageRef &irPos, CVD::BasicImage<CVD::byte> &im, int nRange, bool showStats)
{
  if(!im.in_image_with_border(irPos,nRange))
  {
      return false;
  }

  ImageRef irBest;
  int nBestSSD = 1e6;
  ImageRef square = ImageRef(nRange, nRange);
  ImageRef irBBoxTL = irPos - square;
  ImageRef irBBoxBR = irPos + square;

  for(int x=irBBoxTL.x; x<irBBoxBR.x;x++){
    for(int y=irBBoxTL.y; y<irBBoxBR.y; y++){
      ImageRef irTarget(x,y);
      int nSSD = SSDAtPoint(im, irTarget);
      if(showStats)
        cout << "SSD for point " << irTarget << " : " << nSSD << endl;
      if(nSSD == nBestSSD)
      {
        irPos = ImageRef(-1, -1);
        return false;
      }
      if(nSSD < nBestSSD)
      {
        irBest = irTarget;
        nBestSSD = nSSD;
      }
    }
  }

  irPos = irBest; // always set it, just in case

  if(nBestSSD < mnMaxSSD)
    return true;
  else
    return false;
}

// Define the patch from an input image
void MiniPatch::SampleFromImage(ImageRef irPos, BasicImage<CVD::byte> &im)
{
  assert(im.in_image_with_border(irPos, mnHalfPatchSize));
  CVD::ImageRef irPatchSize( 2 * mnHalfPatchSize + 1 , 2 * mnHalfPatchSize + 1);
  mimOrigPatch.resize(irPatchSize);
  copy(im, mimOrigPatch, mimOrigPatch.size(), irPos - mimOrigPatch.size() / 2);
}

// Static members
int MiniPatch::mnHalfPatchSize = 4;
int MiniPatch::mnRange = 10;
int MiniPatch::mnMaxSSD = GV2.GetInt("MaxSSD", 19999, SILENT);
