/**
 * Date 19/11/12
 * Author Sebastian Bronte <sebastian.bronte@depeca.uah.es>
 * This class inherits MapPoint class to implement specific functionality 
 * related to Non-Rigid Point Handling
 */

#ifndef __NRMAPPOINT_H_
#define __NRMAPPOINT_H_

#include "MapPoint.h"

struct NRMapPoint : public MapPoint
{
  //array of deformable versions for the same point
  std::vector<Vector<3> > v3DefVersions;
};

#endif
