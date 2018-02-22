// -*- c++ *--
// Sebastian Bronte <sebastian.bronte@depeca.uah.es>
// Date 8/12/2011
//
// Interpolation.h
// Declares the Interpolation class
// 
// This is a very simple class to provide a interpolation service;
// this can be replaced with whatever form of interpolation needed. 
// It should take the base points as reference and interpolate those 
// points to the detected ones and vice-versa. Wrap and UnWrap 
// functions provide the interface for the rest of the inherited 
// classes.
//

#ifndef _INTERP_
#define _INTERP_

#include <cvd/byte.h>
#include <TooN/TooN.h>
#include <TooN/wls.h>
#include <TooN/SVD.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "DBases.h"

class NRTrackerData;
class ATANCamera;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_2<K>        Surface;
typedef Surface::Point                           SPoint;
typedef Surface::Face_handle                     Face_handle;
typedef Surface::Vertex_handle                   Vertex_handle;
typedef Surface::Edge_iterator                   Edge_iterator;
typedef Surface::Vertex_iterator                 Vertex_iterator;
typedef Surface::Segment                         Segment;
typedef Surface::Triangle                        Triangle;
typedef Surface::Face_circulator                 Face_circulator;
typedef Surface::Edge_circulator                 Edge_circulator;

class Interpolation
{
  public:
    //class constructor prototype
    Interpolation(ATANCamera & camera, DBases &bases, double * defCoeffs,
                  double * rigidCoeffs, SE3<>& trackerPose);
    virtual ~Interpolation(){}
    // function prototype to apply the interpolation from a deformation
    // base to the detected points
    void Wrap(std::vector<NRTrackerData*> &vTD);
    // function prototype to apply the inverse transformation: from the
    // detected points to the defromation base points
    void UnWrap(std::vector<NRTrackerData*> &vTD);

    //This function is in charge of drawing the current triangulation from the
    //deformation bases
    void drawTriangulation();

    //This function computes the 3D shape from the set of coefficients.
    void compute3DShape();

    //This function projects into the 2D shape structure the pre-computed 3D one...
    void compute2DShape();

    std::vector<unsigned int> & getNeighboursFromPointIndex(unsigned int idx){
      return idx<Vertex_neighbourIdx.size() ? Vertex_neighbourIdx[idx] : emptyVector;
    }

    Vector<3,unsigned int> & getSelectedNeighboursIdx(unsigned int idx){
      return selectedNeighbours[idx];
    }

    Vector<3> & getConstraintsNeighbours(unsigned int idx){
      return constraintsNeighbours[idx];
    }

    std::vector<Vector<2> > m2Shape;
    std::vector<Vector<3> > m3Shape;

    inline int nFaces(){return S.number_of_faces();}

    inline std::pair<Vector<3>, double> & getTemplateDiff(unsigned int idx1,
                                                          unsigned int idx2){
      unsigned int id1 = idx1, id2 = idx2;
      if(idx1>idx2){
        id1 = idx2;
        id2 = idx1;
      }

      std::pair<unsigned int, unsigned int> p(id1,id2);
      if(pairDist.find(p) != pairDist.end()){
        return pairDist[p];
      }else{
        return emptyPair;
      }
    }

  protected:
    ATANCamera & mCamera;
    DBases & mBases; //Reference to the base object

    SE3<> & mse3TrackerPose; //reference to the current pose from the tracker

    Surface S;

    std::map<Surface::Point, std::size_t > Vertex_map;
    std::vector<std::vector<unsigned int> >Vertex_neighbourIdx;
    std::vector<unsigned int> emptyVector;
    std::pair<Vector<3>,double> emptyPair;

    std::vector<Vector<3> > m3ShapeTempl;
    Matrix<> *DistTempl;
    std::map<std::pair<unsigned int,unsigned int>, std::pair<Vector<3>,double> > pairDist;

    double * mdDefCoeffs;
    double * mdRigidCoeffs;

    std::vector<Vector<3, unsigned int> > selectedNeighbours;
    std::vector<Vector<3> > constraintsNeighbours;

    // The following function shouldn't be re-implemented since it is common 
    // for all the inheritance. They should be called in the class constructor, 
    // although some track information must be provided after, so it might be
    // better calling them from the Wrap and Unwrap functions
    // Faces are estimated from Dbases shape points to further estimate the
    // centroids for each detected point to each point in the base.
    // In the case that the faces are extended, the function should be called
    // again (from the MapMaker thread)
    void estimateFacesFromRigidShape();

    // (Only valid for regular meshes) This function estimates the distances to
    // the nearest neighbour point within the template shape. It is supposed
    // that the points on the templates are aligned to the rest of the map,
    // at least on the first frame.
    virtual void estimateTemplateDistances()=0;

    //here the per-point system solution must be implemented
    virtual void centroidFromFaces(std::vector<NRTrackerData*> &vTD)=0;
    virtual void interpolate(NRTrackerData* TD)=0;
    virtual void deInterpolate(NRTrackerData* TD)=0;
};

// Nearest neighbour interpolation, to speed up DIRECT correspondences, 
// no interpolation is applied in principle

class Interpolation_NN: public Interpolation 
{
  public:
    Interpolation_NN(ATANCamera & camera, DBases &bases, double * defCoeffs,
                     double * rigidCoeffs, SE3<>& trackerPose):
      Interpolation(camera, bases, defCoeffs, rigidCoeffs, trackerPose){}
    virtual ~Interpolation_NN(){}
  protected:
    virtual void centroidFromFaces(std::vector<NRTrackerData*> &vTD);
    virtual void interpolate(NRTrackerData* TD);
    virtual void deInterpolate(NRTrackerData* TD);
    virtual void estimateTemplateDistances();
};

// Linear interpolation class

class Interpolation_LN: public Interpolation
{
  public:
    Interpolation_LN(ATANCamera & camera, DBases &bases, double * defCoeffs,
                     double * rigidCoeffs, SE3<>& trackerPose):
      Interpolation(camera, bases, defCoeffs, rigidCoeffs, trackerPose){}
    virtual ~Interpolation_LN(){}
  protected:
    virtual void centroidFromFaces(std::vector<NRTrackerData*> &vTD);
    virtual void interpolate(NRTrackerData* TD);
    virtual void deInterpolate(NRTrackerData* TD);
    virtual void estimateTemplateDistances();

    void checkPointLinearityAndRetrieveNeighbourIdx(std::vector<unsigned int> &vIdx,
                                                    Vector<3, unsigned int> &idx);
    void computeShapeConstraints(Vector<3> &point, std::vector<unsigned int> &vIdx,
                                 Vector<3, unsigned int> &idx, Vector<3> &constraints);
};

// Other interpolation algorithms can be tried out, such as the ones based on
// splines, but I had no time to do so.

#endif
