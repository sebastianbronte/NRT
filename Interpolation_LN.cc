/*
* Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
* Date: 12/12/2011
* Interpolation_LN class definition
*/

#include "Interpolation.h"
#include "NRTrackerData.h"

#include <TooN/wls.h>

using namespace std;

/**
 * computes the centroids from the current set of faces
 * @param vTD
 */
void Interpolation_LN::centroidFromFaces(vector<NRTrackerData*> &vTD)
{
  //look for the facet of every point
  //Then, save in the TrackerData structure the centroid
  unsigned int nbases = mBases.getNBasis();
  for(unsigned int i=0;i<vTD.size();i++)
  {
    NRTrackerData * TD = vTD[i];
    if(TD==NULL)
      break;
    //if already computed, skip this point
    if(TD->isAttached() || !TD->bFound || (TD->bFound && TD->nSearchLevel>1))
      continue;

    Vector<2> p = TD->v2Image;
    SPoint p2 = SPoint(p[0],p[1]);
    Face_handle f = S.locate(p2);
    if(f==NULL || S.is_infinite(f))
    {
      TD->attached = false;
      continue;
    }
    TD->face = f;

    Triangle t = S.triangle(f);
    SPoint v0 = t.vertex(0);
    SPoint v1 = t.vertex(1);
    SPoint v2 = t.vertex(2);
    size_t idx0 = Vertex_map[v0];
    size_t idx1 = Vertex_map[v1];
    size_t idx2 = Vertex_map[v2];

    // Uncomment if you want to try this way of computing the centroids,
    // but I guess it is a bit slower than just taking the data.
    // On the other hand, with this, there is no need to be inside the triangle
    /*Matrix<2,3> A;
    A[0][0]=v0[0];A[1][0]=v0[1];
    A[0][1]=v1[0];A[1][1]=v1[1];
    A[0][2]=v2[0];A[1][2]=v2[1];

    SVD<> wls(A);
    Vector<3> centroid = wls.backsub(makeVector(p[0],p[1]));*/
    Vector<2> vv0, vv1, vv2;
    vv0[0] = v0[0];vv0[1] = v0[1];
    vv1[0] = v1[0];vv1[1] = v1[1];
    vv2[0] = v2[0];vv2[1] = v2[1];
    Vector<2> a0 = vv1 - vv0;
    Vector<2> a1 = vv2 - vv0;
    Vector<2> a2 = p - vv0;
    float d00 = a0*a0;
    float d01 = a0*a1;
    float d11 = a1*a1;
    float d20 = a2*a0;
    float d21 = a2*a1;
    float invDenom = 1.0 / (d00*d11 - d01*d01);
    Vector<3> centroid;
    centroid[1] = (d11*d20-d01*d21)*invDenom;
    centroid[2] = (d00*d21-d01*d20)*invDenom;
    centroid[0] = 1.0 - centroid[1] - centroid[2];

    TD->v3Centroid = centroid;
    TD->attached=true;
    TD->setIndexes(idx0, idx1, idx2);
    TD->vv3BasesPoint.clear();
    //update rigid component
    TD->v3Rigid = TD->v3Centroid[0]*(*mBases.mmRigidShape)[TD->idxs[0]].slice<0,3>()+
                  TD->v3Centroid[1]*(*mBases.mmRigidShape)[TD->idxs[1]].slice<0,3>()+
                  TD->v3Centroid[2]*(*mBases.mmRigidShape)[TD->idxs[2]].slice<0,3>();

    for(unsigned int k=0;k<nbases;k++)
    {
      Vector<3> b0 = mBases.getBasisPoint(k,idx0);
      Vector<3> b1 = mBases.getBasisPoint(k,idx1);
      Vector<3> b2 = mBases.getBasisPoint(k,idx2);
      Vector<3> tmp =  centroid[0]*b0 + centroid[1]*b1 + centroid[2]*b2;
      TD->vv3BasesPoint.push_back(tmp);
    }
    TD->updateDef3D = true;    
  }
}

void Interpolation_LN::estimateTemplateDistances(){
  //m3Shape matrix is supposed to have nPoints x 3.
  //we also supposed the indices aligned  ...
  Vertex_neighbourIdx.clear();
  for(unsigned int i=0;i<m3Shape.size();i++)
  {
    Vector<2> v2 = m2Shape[i];
    SPoint p2(v2[0],v2[1]);
    Vertex_handle vh1 = S.nearest_vertex(p2);
    Edge_circulator ec1 = S.incident_edges(vh1),done(ec1);
    Segment s;
    vector<unsigned int> vIdx = vector<unsigned int>();
    do{    
      if (ec1 != 0 && ! S.is_infinite(ec1))
      { 
        s = S.segment(ec1);
        SPoint p = s.source();

        unsigned int j = Vertex_map[p];
        vIdx.push_back(j);
      }else{
        break;
      }
    }while(++ec1!=done);
    Vertex_neighbourIdx.push_back(vIdx);
  }

  selectedNeighbours.clear();
  constraintsNeighbours.clear();

  for(unsigned int i=0;i<m3Shape.size();i++)
  {
    vector<unsigned int> &vIdx = Vertex_neighbourIdx[i];
    Vector<3,unsigned int> idx = Zeros;
    checkPointLinearityAndRetrieveNeighbourIdx(vIdx,idx);
    // be carefull with all the indices set to zero
    selectedNeighbours.push_back(idx);
    //after the execution of this function, we have all the points ready and
    //their linearity checked, so we have the correct neighbourhood ready
    //to compute the constraints
    Vector<3> constraints = Zeros;
    computeShapeConstraints(m3Shape[i],vIdx,idx,constraints);
    constraintsNeighbours.push_back(constraints);
  }
}

void Interpolation_LN::checkPointLinearityAndRetrieveNeighbourIdx(vector<unsigned int> &vIdx, Vector<3, unsigned int> &idx){
  static const double epsilon = 0.001;
  Vector<3> v[3], vtmp, vtmp2;
  unsigned int idxtmp, idxtmp2;
  double det=0;
  idx = Zeros;

  if (vIdx.size() < 3)
    return;

  for (unsigned int k = 0; k<vIdx.size() && abs(det) < epsilon;k++){
    if (k < 3){
      idx[k] = k;
      v[k] = m3Shape[vIdx[k]];
      if (k==2)
        det = v[0][0]*v[1][1]*v[2][2]+v[0][2]*v[1][0]*v[2][1]+v[0][1]*v[1][2]*v[2][0]-v[0][2]*v[1][1]*v[2][0]-v[0][0]*v[1][2]*v[2][1]-v[0][1]*v[1][0]*v[2][2];
    }else{
      int trial = 0;
      vtmp2 = v[trial];
      idxtmp2 = idx[trial];
      vtmp = m3Shape[vIdx[k]];
      idxtmp = k;
      do{
        vtmp2 = v[trial];
        idxtmp2 = idx[trial];
        v[trial] = vtmp;
        idx[trial] = idxtmp;
        //Krammer rule, for 3x3 matrix, faster than Cholesky decomposition?
        det = v[0][0]*v[1][1]*v[2][2]+v[0][2]*v[1][0]*v[2][1]+v[0][1]*v[1][2]*v[2][0]-v[0][2]*v[1][1]*v[2][0]-v[0][0]*v[1][2]*v[2][1]-v[0][1]*v[1][0]*v[2][2];
        if (abs(det) < epsilon){
          v[trial] = vtmp2;
          idx[trial] = idxtmp2;
        }
        trial++;
      }while(abs(det) < epsilon && trial < 3);       
    }
  }

  // non linearly dependent
  if(abs(det) < epsilon)
    idx = Zeros;
}

void Interpolation_LN::computeShapeConstraints(Vector<3> &point, vector<unsigned int> &vIdx, Vector<3, unsigned int> &idx, Vector<3> &constraints)
{ 
  if (vIdx.size() < 3 || (idx[0] == 0 && idx[1] == 0 && idx[2] == 0))
    return;

  Matrix<3> base;
  Vector<3> v0=m3Shape[vIdx[idx[0]]];
  Vector<3> v1=m3Shape[vIdx[idx[1]]];
  Vector<3> v2=m3Shape[vIdx[idx[2]]];
  base[0][0]=v0[0];base[0][1]=v1[0];base[0][2]=v2[0];
  base[1][0]=v0[1];base[1][1]=v1[1];base[1][2]=v2[1];
  base[2][0]=v0[2];base[2][1]=v1[2];base[2][2]=v2[2];
  SVD<> svd(base);
  constraints = svd.backsub(point);
}


/*
 * interpolate
 * @brief This function is in charge of providing correspondences
 *        between detected and deformation base points using linear interpolation.
 *
 * @param p1 original point (detected)
 * @param p2 destination points (basis points)
 * @param centroids centroids of p1 in face
 * @param f face in which p1 is located
 *
 */
void Interpolation_LN::interpolate(NRTrackerData* TD)
{
  //Nothing is done on purpose
}

/*
 * deInterpolate
 * @brief This function is in charge of provide correspondences
 *        between detected and deformation base points using linear interpolation.
 *
 * @param p1 original points (basis points)
 * @param p2 destination point (detected)
 * @param centroids centroids of p1 in face
 * @param f face in which p1 is located
 *
 */
void Interpolation_LN::deInterpolate(NRTrackerData* TD)
{
  Vector<3> centroid = TD->v3Centroid;
  if(!TD->isAttached()){
      return;
  }

  size_t idx0 = TD->idxs[0]; //Original base index of point v0
  size_t idx1 = TD->idxs[1]; //Original base index of point v1
  size_t idx2 = TD->idxs[2]; //Original base index of point v2

  Vector<2> p2D = Zeros;
  p2D += centroid[0]*m2Shape[idx0];
  p2D += centroid[1]*m2Shape[idx1];
  p2D += centroid[2]*m2Shape[idx2];
  TD->v2Image = p2D;
}
