// -*- c++ -*-
// Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
// Date: 5 May 2011
// NRTrackerData.h
//
// Tracker Data handler for Non Rigid SfM tracker class. Inherited from
// TrackerData original handler, but adapted to handle deformation coefficients,
// bases and implement the rest of equations required
//

#ifndef __NRTRACKERDATA__
#define __NRTRACKERDATA__

#include "TrackerData.h"
#include "DBases.h"
#include "ATANCamera.h"
#include "NRMapPoint.h"

#include <sstream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <opencv/cv.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>        Surface;
typedef Surface::Face_handle                     Face_handle;

class DBases;

struct NRTrackerData: public TrackerData
{
  NRTrackerData(NRMapPoint *pMapPoint, DBases * base = NULL,
               const int pointidx = 0, double *coefs = NULL):
               TrackerData(pMapPoint),vDefCoefs(coefs){
    //all the data can be dimensioned since we have the bases already given.
    Base = base;
    unsigned int nbases = Base==NULL?0:Base->getNBasis();

    m2CoefDerivs = nbases==0?NULL:new double[nbases*2];
    m2kJacobian = nbases==0?NULL:new double[nbases*2];
    mp26kJacobian = nbases==0?NULL:new double[2*(nbases+6)];

    //old version index
    pidx = pointidx;

    //new version indexes. to be initialized later
    idxs = -1*Ones;

    v3Def = Zeros;
    v3Rigid = Zeros;
    v3Map = Zeros;
    v3MapPrev = Zeros;

    updateDef3D=nbases==0?false:true;

    v3Centroid = Zeros;
    attached = false;
    face=NULL;

    oldTracks = std::vector<Vector<2> >();

    mDescriptor = NULL;
  } //constructor

  //Destructor. Some members are dynamicaly built, then they need to be deleted
  ~NRTrackerData()
  {
    delete [] m2CoefDerivs;
    delete [] m2kJacobian;
    delete [] mp26kJacobian;
    if(mDescriptor) delete [] mDescriptor;
  } 

  // Deformation coeficient derivs. Dimensionality dynamically constructed
  // from the number of basis on the base
  double *m2CoefDerivs;
  // Deformation coeficients to apply to each model (map) point. Dynamically
  // constructed since the length of the coeficients is undefined as well
  double *vDefCoefs; 

  Vector<3> v3Def; //deformation component of the 3D point
  Vector<3> v3Rigid; //rigid component of the 3D point (average shape)

  Vector<3> v3Map;
  Vector<3> v3MapPrev;

  Vector<3> v3Centroid;
  std::vector<Vector<3> > vv3BasesPoint;
  Face_handle face;
  bool attached;

  bool isAttached(){return attached;}

  DBases * Base; //reference to the loaded base in the program.

  //old version index
  unsigned int pidx; //index of the point in the bases

  //new version indexes
  Vector<3,int> idxs;
  std::vector<unsigned int> vNeighbourIdx;

  // The part of the Jacobian responsible to measure the error associated
  // to the coeficient estimation
  double *m2kJacobian; 

  double *mp26kJacobian;
  bool updateDef3D;

  Vector<6> priorJacobian;

  std::vector<Vector<2> > oldTracks;

  unsigned char * mDescriptor;
  cv::KeyPoint mKPoint;

  inline void Deform3DTerm()
  {
    v3Def = Zeros;
    unsigned int nbases = Base==NULL?0:Base->getNBasis();
    if(nbases)
    {
      if(idxs[0]==-1 || idxs[1]==-1 || idxs[2]==-1)
      {
        if(pidx > Base->getNPoints())
          return;
        for (unsigned int i = 0; i < nbases; i++)
        {
          //get one of the basis->take a slice from the whole base.
          Vector <3> basis_i = Base->getBasisPoint(i,pidx);
          //calculate the part due to deformation for each 3D point
          v3Def += vDefCoefs[i]*basis_i;
        }
      }else{
        for (unsigned int i = 0; i < nbases; i++)
          v3Def+= vDefCoefs[i]*vv3BasesPoint[i];
      }
      //std::cout << "3d point def " << v3Def << " 3d point rigid " << v3Rigid
      //          << " final point " << v3Rigid + v3Def << std::endl;
    }

    updateDef3D=false;
  }

  // Projection of a non rigid
  inline void Project(const SE3<> &se3CFromW, ATANCamera &Cam)
  {
    bInImage = bPotentiallyVisible = false;
    //temporal hack: do not compute v3Def if it's already computed.
    //It can be clearly seen that if the coeficients don't get updated,
    //compute them every time, it's a huge waste of time
    if(updateDef3D)
      Deform3DTerm();
    if(Base!=NULL)
    {
      //the final point is formed by the rigid and deformable part
      v3Map = v3Def + v3Rigid;
    }
    Point.v3WorldPos = v3Map;
    v3Cam = se3CFromW * Point.v3WorldPos;
    if(v3Cam[2] < 0.001)
      return;
    //std::cout << "v3Def " << v3Def << " v3Rigid " << v3Rigid << " total "
    //          << Point.v3WorldPos << " transformed " << v3Cam << std::endl;

    v2ImPlane = project(v3Cam);

    if(v2ImPlane*v2ImPlane > Cam.LargestRadiusInImage() * Cam.LargestRadiusInImage())
      return;
    v2Image = Cam.Project(v2ImPlane);
    //std::cout << " v2Image " << v2Image << " v2Found " << v2Found 
    //          << " err " << v2Image - v2Found << std::endl;

    if(Cam.Invalid())
      return;

    if(v2Image[0] < 0 || v2Image[1] < 0 || v2Image[0] > irImageSize[0] || v2Image[1] > irImageSize[1])
      return;
    bInImage = true;
  }

  // Get the projection derivatives (depend only on the camera.)
  // This is called Unsafe because it depends on the camera caching
  // results from the previous projection:
  // Only do this right after the same point has been projected!
  inline void GetDerivsUnsafe(ATANCamera &Cam)
  {
    m2CamDerivs = Cam.GetProjectionDerivs();
  }

  // Jacobian of projection W.R.T. the camera position
  // I.e. if  p_cam = SE3Old * p_world,
  //         SE3New = SE3Motion * SE3Old
  
  inline void CalcJacobian()
  {
    if(mp26kJacobian!=NULL)
    {
      // If more terms are required in this jacobian, just add them with
      // 6+Base->getNBasis()
      Matrix <Dynamic,Dynamic,double,Reference::RowMajor> m26kJacobian(mp26kJacobian,2,6);
      TrackerData::CalcJacobian();
      m26kJacobian.slice<0,0,2,6>() = m26Jacobian;
    }
    else
    {
      TrackerData::CalcJacobian();
    }
    // In case the derivatives of the coefficients are considered,
    // the 6+NBasis positions in m26kJacobian should be filled in.
  }

  // This function had changed its prototype since we need the pose
  // to estimate the derivative of each one of the basis.
  inline void CalcJacobian(const SE3<> & se3CFromW, const SE3<> & se3CFromWPrev)
  {
    this->CalcJacobian();
    //compute the difference between previous and current poses
    SE3<> poseDiff = se3CFromW * se3CFromWPrev.inverse();
    priorJacobian = 2*SE3<>::ln(poseDiff);
  }

  // Sometimes, within tracker class, instead of a full reprojection,
  // the position is updated linearly.
  inline void LinearUpdate(const Vector<> &v6k)
  {
    v2Image += m26Jacobian * v6k.slice<0,6>();
  }

  inline void ProjectAndDerivs(SE3<> &se3, ATANCamera &Cam)
  {
    Project(se3, Cam);
    if(bFound)
      GetDerivsUnsafe(Cam);
  }

  inline void setIndexes(int idx1, int idx2, int idx3)
  {
    idxs = makeVector(idx1, idx2, idx3);
  }

};

#endif
