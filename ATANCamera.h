// *-* c++ *-*
// Copyright 2008 Isis Innovation Limited

// N-th implementation of a camera model
// GK 2007
// Evolved a half dozen times from the CVD-like model I was given by
// TWD in 2000
// 
// This one uses the ``FOV'' distortion model of
// Deverneay and Faugeras, Straight lines have to be straight, 2001
//
// BEWARE: This camera model caches intermediate results in member variables
// Some functions therefore depend on being called in order: i.e.
// GetProjectionDerivs() uses data stored from the last Project() or UnProject()
// THIS MEANS YOU MUST BE CAREFUL WITH MULTIPLE THREADS
// Best bet is to give each thread its own version of the camera!
//
// Camera parameters are stored in a GVar, but changing the gvar has no effect
// until the next call to RefreshParams() or SetImageSize().
//
// Pixel conventions are as follows:
// For Project() and Unproject(),
// round pixel values - i.e. (0.0, 0.0) - refer to pixel centers
// I.e. the top left pixel in the image covers is centered on (0,0)
// and covers the area (-.5, -.5) to (.5, .5)
//
// Be aware that this is not the same as what opengl uses but makes sense
// for acessing pixels using ImageRef, especially ir_rounded.
//
// What is the UFB?
// This is for projecting the visible image area
// to a unit square coordinate system, with the top-left at 0,0,
// and the bottom-right at 1,1
// This is useful for rendering into textures! The top-left pixel is NOT
// centered at 0,0, rather the top-left corner of the top-left pixel is at 
// 0,0!!! This is the way OpenGL thinks of pixel coords.
// There's the Linear and the Distorting version - 
// For the linear version, can use 
// glMatrixMode(GL_PROJECTION); glLoadIdentity();
// glMultMatrix(Camera.MakeUFBLinearFrustumMatrix(near,far));
// To render un-distorted geometry with full frame coverage.
//

// A couple of more functions are added to give the deformation constants
// w.r.t the camera, or directly undistording pixel coordinates.

#ifndef __ATAN_CAMERA_H
#define __ATAN_CAMERA_H

#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <cmath>
using namespace TooN;
#include <cvd/vector_image_ref.h>
#include <gvars3/gvars3.h>

#include "opencv2/core/core.hpp"

#define NUMTRACKERCAMPARAMETERS 5

class CameraCalibrator;
class CalibImage;

// The parameters are:
// 0 - normalized x focal length
// 1 - normalized y focal length
// 2 - normalized x offset
// 3 - normalized y offset
// 4 - w (distortion parameter)

class ATANCamera {
 public:
  ATANCamera(std::string sName);
  ATANCamera(std::string sName, CVD::ImageRef irSize, Vector<NUMTRACKERCAMPARAMETERS> vParams);

  // Image size get/set: updates the internal projection params to that target image size.
  void SetImageSize(Vector<2> v2ImageSize);
  inline void SetImageSize(CVD::ImageRef irImageSize) {SetImageSize(vec(irImageSize));}
  inline Vector<2> GetImageSize() {return mvImageSize;}
  void RefreshParams();
  
  // Various projection functions
  Vector<2> Project(const Vector<2>& camframe); // Projects from camera z=1 plane to pixel coordinates, with radial distortion
  inline Vector<2> Project(CVD::ImageRef ir) { return Project(vec(ir)); }
  Vector<2> UnProject(const Vector<2>& imframe); // Inverse operation
  inline Vector<2> UnProject(CVD::ImageRef ir)  { return UnProject(vec(ir)); }

  Vector<2> UFBProject(const Vector<2>& camframe);
  Vector<2> UFBUnProject(const Vector<2>& camframe);
  inline Vector<2> UFBLinearProject(const Vector<2>& camframe);
  inline Vector<2> UFBLinearUnProject(const Vector<2>& fbframe);
  
  Matrix<2,2> GetProjectionDerivs(); // Projection jacobian
  
  inline bool Invalid() {  return mbInvalid;}
  inline double LargestRadiusInImage() {  return mdLargestRadius; }
  inline double OnePixelDist() { return mdOnePixelDist; }
  
  // OpenGL helper function
  Matrix<4> MakeUFBLinearFrustumMatrix(double near, double far);

  // Feedback for Camera Calibrator
  double PixelAspectRatio() { return mvFocal[1] / mvFocal[0];}
  
  
  // Useful for gvar-related reasons (in case some external func tries to read the camera params gvar, and needs some defaults.)
  static const Vector<NUMTRACKERCAMPARAMETERS> mvDefaultParams;
  
  inline Vector<2> ComputeCoefDerivs(const Vector<3>& Base, const Vector<3>& v3Point, const SE3<> &se3CamFromWorld);
  inline Vector<2> GetBasisTermForCoeffs(const Vector<3>& Base, const Vector<2>& v2Point, const Matrix<3>& R); //to compute term for matrix A in order to get the coeffs
  inline Vector<2> GetRestTermForCoeffs(const Vector<3>& T, const Vector<2>& v2Point); //to compute term for matrix B in order to get the coeffs
  inline Vector<2> GetRestTermForCoeffs(const Vector<3>& T, const Vector<2>& v2Point, const Vector<3>& v3RigidPoint, const Matrix<3> & R);
  inline Vector<2> Distort(const Vector<2>& vIm); //distord one point in camera plane and gives the result in the same plane
  inline Vector<2> UnDistort(const Vector<2>& vIm); //undistord one point in the camera plane and gives the result in the same plane

  inline cv::Mat & getocvCameraMatrix(){ return mocvCameraMatrix; }
  
 protected:
  GVars3::gvar3<Vector<NUMTRACKERCAMPARAMETERS> > mgvvCameraParams; // The actual camera parameters
  
  Matrix<2, NUMTRACKERCAMPARAMETERS> GetCameraParameterDerivs();
  void UpdateParams(Vector<NUMTRACKERCAMPARAMETERS> vUpdate);
  void DisableRadialDistortion();
  
  // Cached from the last project/unproject:
  Vector<2> mvLastCam;      // Last z=1 coord
  Vector<2> mvLastIm;       // Last image/UFB coord
  Vector<2> mvLastDistCam;  // Last distorted z=1 coord
  double mdLastR;           // Last z=1 radius
  double mdLastDistR;       // Last z=1 distorted radius
  double mdLastFactor;      // Last ratio of z=1 radii
  bool mbInvalid;           // Was the last projection invalid?
  
  // Cached from last RefreshParams:
  double mdLargestRadius; // Largest R in the image
  double mdMaxR;          // Largest R for which we consider projection valid
  double mdOnePixelDist;  // z=1 distance covered by a single pixel offset (a rough estimate!)
  double md2Tan;          // distortion model coeff
  double mdOneOver2Tan;   // distortion model coeff
  double mdW;             // distortion model coeff
  double mdWinv;          // distortion model coeff
  double mdDistortionEnabled; // One or zero depending on if distortion is on or off.
  Vector<2> mvCenter;     // Pixel projection center
  Vector<2> mvFocal;      // Pixel focal length
  Vector<2> mvInvFocal;   // Inverse pixel focal length
  Vector<2> mvImageSize;  
  Vector<2> mvUFBLinearFocal;
  Vector<2> mvUFBLinearInvFocal;
  Vector<2> mvUFBLinearCenter;
  Vector<2> mvImplaneTL;   
  Vector<2> mvImplaneBR;
  
  cv::Mat mocvCameraMatrix;
  
  void computeocvCameraMatrix();
  
  // Radial distortion transformation factor: returns ration of distorted / undistorted radius.
  inline double rtrans_factor(double r)
  {
    if(r < 0.001 || mdW == 0.0)
      return 1.0;
    else 
      return (mdWinv* atan(r * md2Tan) / r);
  }

  // Inverse radial distortion: returns un-distorted radius from distorted.
  inline double invrtrans(double r)
  {
    if(mdW == 0.0)
      return r;
    return(tan(r * mdW) * mdOneOver2Tan);
  }
  
  std::string msName;

  // friend declarations allow access to calibration jacobian and camera update function.
  friend class CameraCalibrator;
  friend class CalibImage;
};

// Some inline projection functions:
inline Vector<2> ATANCamera::UFBLinearProject(const Vector<2>& camframe)
{
  Vector<2> v2Res;
  v2Res[0] = camframe[0] * mvUFBLinearFocal[0] + mvUFBLinearCenter[0];
  v2Res[1] = camframe[1] * mvUFBLinearFocal[1] + mvUFBLinearCenter[1];
  return v2Res;
}

inline Vector<2> ATANCamera::UFBLinearUnProject(const Vector<2>& fbframe)
{
  Vector<2> v2Res;
  v2Res[0] = (fbframe[0] - mvUFBLinearCenter[0]) * mvUFBLinearInvFocal[0];
  v2Res[1] = (fbframe[1] - mvUFBLinearCenter[1]) * mvUFBLinearInvFocal[1];
  return v2Res;
}

inline Vector<2> ATANCamera::ComputeCoefDerivs(const Vector<3>& Base,
                                           const Vector<3>& v3Point,
                                           const SE3<> &se3CamFromWorld)
{
  Vector<2> v2Res;
  Matrix<3> R = se3CamFromWorld.get_rotation().get_matrix();
  Vector<3> T = se3CamFromWorld.get_translation();

  double factor = rtrans_factor(mdLastR);
  double rxX = R[0]*v3Point;
  double rxB = R[0]*Base;
  double ryX = R[1]*v3Point;
  double ryB = R[1]*Base;
  double rzX = R[2]*v3Point;
  double rzB = R[2]*Base;

  double lambda = rzX+T[2];

  v2Res[0] = (mvFocal[0]*factor*rxB)/lambda-(rzB*mvFocal[0]*factor*(rxX+T[0]))/(lambda*lambda);
  v2Res[1] = (mvFocal[1]*factor*ryB)/lambda-(rzB*mvFocal[1]*factor*(ryX+T[1]))/(lambda*lambda);

  return v2Res;
}

inline Vector<2> ATANCamera::GetBasisTermForCoeffs(const Vector<3>& Base, const Vector<2>& v2Point, const Matrix<3>& R)
{
  Vector<2> res;
  res[0] = ((v2Point[0]-mvCenter[0])*R[2]-mvFocal[0]*R[0])*Base;
  res[1] = ((v2Point[1]-mvCenter[1])*R[2]-mvFocal[1]*R[1])*Base;
  return res;
}

inline Vector<2> ATANCamera::GetRestTermForCoeffs(const Vector<3>& T, const Vector<2>& v2Point)
{
  Vector<2> res;
  res[0] = mvFocal[0]*T[0]-T[2]*(v2Point[0]-mvCenter[0]);
  res[1] = mvFocal[1]*T[1]-T[2]*(v2Point[1]-mvCenter[1]);
  return res;
}

inline Vector<2> ATANCamera::GetRestTermForCoeffs(const Vector<3>& T, const Vector<2>& v2Point, const Vector<3>& v3RigidPoint, const Matrix<3> & R)
{
  Vector<3> rotRigidPoint = R*v3RigidPoint;
  Vector<2> res;
  res[0] = mvFocal[0]*(T[0]+rotRigidPoint[0])-(T[2]+rotRigidPoint[2])*(v2Point[0]-mvCenter[0]);
  res[1] = mvFocal[1]*(T[1]+rotRigidPoint[1])-(T[2]+rotRigidPoint[2])*(v2Point[1]-mvCenter[1]);
  return res;
}

inline Vector<2> ATANCamera::Distort(const Vector<2>& vIm)
{
  Vector<2> vInter;
  vInter[0] = (vIm[0]-mvCenter[0])*mvInvFocal[0];
  vInter[1] = (vIm[1]-mvCenter[1])*mvInvFocal[1];

  return Project(vInter);
}

inline Vector<2> ATANCamera::UnDistort(const Vector<2>& vIm)
{
  Vector<2> vRes, vInter;
  vInter = UnProject(vIm);
  vRes[0] = vInter[0]*mvFocal[0]+mvCenter[0];
  vRes[1] = vInter[1]*mvFocal[1]+mvCenter[1];
  return vRes;
}

#endif
