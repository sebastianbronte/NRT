/**
 * Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
 * Date: 8/12/2011
 *
 * Interpolation.cpp
 * Class which implements the different interpolation 
 * methods from inherited Interpolation classes.
 */

#include "ATANCamera.h"
#include "Interpolation.h"
#include "NRTrackerData.h"
#include "OpenGL.h"

#include <cvd/vision.h>
#include <algorithm>

#include <sys/time.h>

using namespace CVD;
using namespace std;

/**
 * @brief Interpolation constructor
 * @param camera Reference to the camera for projections or 2D point operations
 * @param bases Reference to the object containing the bases.
 * @param defCoeffs the pointer to the set of coefficients
          (a shallow copy will be done).
 * @param rigidCoeffs the pointer to the set of rigid coefficients
          (a shallow copy will be done).
 * @param trackerPose reference to the pose contained by the tracker
 */
Interpolation::Interpolation(ATANCamera & camera, DBases &bases,
        double * defCoeffs, double * rigidCoeffs, SE3<>& trackerPose)
       :m2Shape(vector<Vector<2> >(bases.getNPoints())),
        m3Shape(vector<Vector<3> >(bases.getNPoints())),
        mCamera(camera),mBases(bases),mse3TrackerPose(trackerPose),
        mdDefCoeffs(defCoeffs), mdRigidCoeffs(rigidCoeffs)
{
  emptyVector.clear();
  emptyPair = pair<Vector<3>,double>(Zeros, 0.0);
}

/**
 * estimateFacesFromBases
 * @brief This function gives an interface from this program to the external
 *        libraries in charge of working out the set of faces from a given
 *        deformation base.
 */
void Interpolation::estimateFacesFromRigidShape()
{
  static bool alreadyComputed = false;
  if(alreadyComputed)
    return;
  //Adapt from Vector to the type required by CGAL
  //do the triangulation from the rigid shape
  for(unsigned int i=0;i<m2Shape.size();i++) {
    //Checking before introducing? ->see CGAL/examples/Triangulation_3/find_conflicts_3.cpp
    Vector<2> b = m2Shape[i];
    SPoint p = SPoint(b[0],b[1]);
    S.insert(p);
    Vertex_map[p] = i;
  }
  assert(S.is_valid());
  estimateTemplateDistances();
  alreadyComputed=true;
}

/**
 * estimateTemplateDistances
 * @brief From the set of points on the template stored in this class, the
 * distances between the closest neighbours points are measured and stored.
 */
void Interpolation::estimateTemplateDistances(){
  //m3ShapeTempl matrix is supposed to have nPoints x 3.
  //we also supposed the indices aligned  ...
  Vertex_neighbourIdx.clear();
  m3ShapeTempl.clear();
  m3ShapeTempl.insert(m3ShapeTempl.begin(),m3Shape.begin(),m3Shape.end());
  for(unsigned int i=0;i<m3Shape.size();i++)
  {
    Vector<3> v3 = m3Shape[i];
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
        unsigned int idx1=i, idx2 = j;
        if (i > j) {
          idx1 = j;
          idx2 = i;
        }

        pair<unsigned int,unsigned int> pidx(idx1, idx2);
        Vector<3> diff = v3 - m3Shape[j];
        double dist = diff*diff;
        pair<Vector<3>,double> diffDist(diff,dist);
        pairDist[pidx]=diffDist;        
      }
    }while(++ec1!=done);
    Vertex_neighbourIdx.push_back(vIdx);
  }
}

/**
 * Wrap
 * @brief This function gives an common interface to interpolate
 *        the points using each class specific interpolation function.
 * @param vTD vector of Tracking points
 */
void Interpolation::Wrap(vector<NRTrackerData*> &vTD)
{
  // I am not very happy checking this all the time.
  static bool alreadyComputed = false;
  if(alreadyComputed)
    return;
  centroidFromFaces(vTD);
  alreadyComputed = true;
}

/**
 * UnWrap
 * @brief This function gives an common interface to de-interpolate
 *        the points using each class specific interpolation function.
 * @param vTD vector of Tracking points
 */
void Interpolation::UnWrap(vector<NRTrackerData*> &vTD)
{
  compute3DShape();
  compute2DShape();
}

/**
 * drawTriangulation
 * @brief Draws the triangulation generated in this class.
 */
void Interpolation::drawTriangulation()
{

  //First draw the lines
  glLineWidth(1);
  glEnable(GL_BLEND);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glBegin(GL_LINES);

  size_t idx1, idx2;
  SPoint v0, v1;
  Vector<2> V0, V1;
  Segment s;
  for(Edge_iterator it = S.edges_begin(); it!=S.edges_end();it++)
  {
    s = S.segment(it);
    v0 = s.vertex(0);
    v1 = s.vertex(1);
    if(s.is_degenerate() || S.is_infinite((*it)))
    {
      continue;
    }
    idx1 = Vertex_map[v0];
    idx2 = Vertex_map[v1];

    glColor3f(0,0.8,0.8);
    V0[0] = m2Shape[idx1][0]; V0[1] = m2Shape[idx1][1];
    glVertex(V0);
    glColor3f(0,0.8,0.8);
    V1[0] = m2Shape[idx2][0]; V1[1] = m2Shape[idx2][1];
    glVertex(V1);
  }
  glEnd();
  glDisable(GL_BLEND);

  //Then draw the points
  glPointSize(2);
  glEnable(GL_BLEND);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glBegin(GL_POINTS);
  for(unsigned int i=0;i<m2Shape.size();i++) {
    glColor3f(0,1,1);
    glVertex(m2Shape[i]);
  }
  glEnd();
  glDisable(GL_BLEND);
}

/**
 * compute3DShape
 * @brief This function computes the 3D shape from the set of coefficients.
 * No pose information is required. It just adds to the rigid shape the
 * non-rigid one by multiplying the set of coeffs to the basis components.
 */
void Interpolation::compute3DShape()
{
  unsigned int nbases = mBases.getNBasis();
  Vector<Dynamic,double,Reference> vDefCoefs(mdDefCoeffs,nbases);
  Vector<Dynamic,double,Reference> vDefCoefs_rigid(mdRigidCoeffs,nbases);

  for(unsigned int p=0;p<mBases.getNPoints();p++)
  {
    Vector<3> pTmp = (*mBases.mmRigidShape)[p].slice<0,3>();
    for(unsigned int k=0;k<mBases.getNBasis();k++)
    {
        Vector<3> b_i = mBases.getBasisPoint(k,p);
        double l_i =mdDefCoeffs[k];
        pTmp += b_i*l_i;
    }

    m3Shape[p] = project(mBases.mse3RelativePose*unproject(pTmp));
  }
}

/**
 * compute2DShape
 * @brief This function projects into the 2D shape structure the pre-computed
 * 3D one.
 */
void Interpolation::compute2DShape()
{
  // Take into account 2 poses: the current one from the tracker, and the one
  // from the map.
  SE3<> finalPose = mse3TrackerPose;

  Vector<3> p3Tmp;
  Vector<2> p2Tmp, p2Proj;

  for(unsigned int p=0; p<m3Shape.size();p++)
  {
    p3Tmp = finalPose*m3Shape[p];

    p2Tmp = project(p3Tmp);

    p2Proj = mCamera.Project(p2Tmp);
    m2Shape[p] = p2Proj;
  }
  estimateFacesFromRigidShape();
}

/**
 * centroidFromFaces
 * @brief This function is in charge of computing the centroids of the detected
 *        point in the context of the Nearest Neighbour algorithm,
 *        i.e. selecting the closest detected point.
 * @param vTD vector of NRTrackerData
 */
void Interpolation_NN::centroidFromFaces(vector<NRTrackerData*> &vTD)
{
  // Take the point
  // Locate the facet
  // Estimate its centroid acording to the NN (Nearest Neighbour)

  // Pix for very sparse points, if closer, this threshold should be lowered.
  double maxdist = 5;

  // Before the pixel by pixel data association, take out the average for each
  // of the data sources, and then, align the data ...

  unsigned int nDet = vTD.size();
  Matrix<> det(nDet,2);
  det = Zeros;
  Vector<2> avg = Zeros;
  for (unsigned int i = 0; i < nDet; i++)
  {
    det[i] = vTD[i]->v2Found;
    avg += vTD[i]->v2Found;
  }

  avg /=nDet;

  for (unsigned int i = 0; i < nDet; i++)
  {
    det[i] -= avg;
  }

  unsigned int nShape = m2Shape.size();
  Matrix<> shape2D2(nShape,2);
  avg = Zeros;
  for (unsigned int i = 0; i < nShape; i++)
  {
    shape2D2[i] = m2Shape[i];
    avg += m2Shape[i];
  }

  avg /= nShape;

  for (unsigned int i = 0; i < nShape; i++){
    shape2D2[i] -= avg;
  }

  //the last ones from the vTD could be possibly corrupted. Avoiding them. 
  for(unsigned int p=0;p<nShape;p++)
  {
    NRTrackerData *TD = vTD[p];

    if(TD==NULL){
      continue;
    }

    TD->attached = false;

    set<int> alreadyUsedIdx;
    alreadyUsedIdx.clear();
    //look for the point on the m2Shape that corresponds to this seen point.
    bool found = false;
    for (unsigned int i = 0; i<m2Shape.size();i++){
     if ((abs(det[p][0]-shape2D2[i][0]) < maxdist) &&
          (abs(det[p][1]-shape2D2[i][1]) < maxdist)){
        if (alreadyUsedIdx.empty() || alreadyUsedIdx.find(i) != alreadyUsedIdx.end()){
          TD->pidx = i;
          alreadyUsedIdx.insert(i);
          found = true;
          break;
        }
      }
    }

    if (!found){
      continue;
    }

    TD->attached=true;

    TD->vv3BasesPoint.clear();
    //update rigid component
    unsigned int idx = TD->pidx;
    TD->v3Rigid = (*mBases.mmRigidShape)[p].slice<0,3>();

    unsigned int nbases = mBases.getNBasis();

    //update non rigid component
    for(unsigned int k=0;k<nbases;k++){
      Vector<3> b = mBases.getBasisPoint(k,idx);
      TD->vv3BasesPoint.push_back(b);
    }

    TD->updateDef3D = true;

    m3Shape[p] = mse3TrackerPose * vTD[p]->v3Map;
    m2Shape[p] = vTD[p]->v2Image;
  }

  estimateFacesFromRigidShape();
}

/**
 * interpolate
 * @brief This function is in charge of providing direct correspondences
 *        between detected and deformation base points. In a Nearest Neighbour
 *        context, the correspondence is inmediate and already filled elsewhere
 * @param TD The pointer to the NRTrackerData structure holding track information.
 */
void Interpolation_NN::interpolate(NRTrackerData* TD)
{
  // Left empty on purpose.
  // There is no need of interpolation in a nearest neighbour.
}

/**
 * deInterpolate
 * @brief This function provides direct correspondences from base points to
 *        detected points. In a Nearest Neighbour implementation, there is no
 *        need to do anything.
 * @param TD The pointer to the NRTrackerData structure holding track information.
 */
void Interpolation_NN::deInterpolate(NRTrackerData* TD)
{
  // Left empty on purpose.
  // There is no need of interpolation in a nearest neighbour.
}

/**
 * estimateTemplateDistances
 * @brief In a Nearest Neighbour context the base implementation of this function is enough.
 */
void Interpolation_NN::estimateTemplateDistances(){
  Interpolation::estimateTemplateDistances();
}

