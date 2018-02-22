// -*- c++ -*-
//Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
//Date: 5 May 2011
//DBase.cc
//
//DBase class definition file. It handles loading deformable bases from disk
//and its operation.
//

#include "DBases.h"

#include <gvars3/instances.h>
#include <gvars3/GStringUtil.h>

#include <sstream>
#include <fstream>

using namespace std;
using namespace GVars3;

DBases::DBases(unsigned int np, unsigned int nb):ncomponents(np),nbasis(nb),
                                                 DBase(NULL),
                                                 mpvBaseComponentEnergy(NULL),
                                                 mpvBaseWeights(NULL),
                                                 normalized(false),
                                                 mbRigidShapeLoadedFromFile(false)
{
  //initialize the relative Map-Bases transformation
  Matrix<3> R;
  R[0] = makeVector(GV3::get<double>("r00b",1),GV3::get<double>("r01b",0),GV3::get<double>("r02b",0));
  R[1] = makeVector(GV3::get<double>("r10b",0),GV3::get<double>("r11b",1),GV3::get<double>("r12b",0));
  R[2] = makeVector(GV3::get<double>("r20b",0),GV3::get<double>("r21b",0),GV3::get<double>("r22b",1));

  Vector<3> T = makeVector(GV3::get<double>("t0b",0),GV3::get<double>("t1b",0),GV3::get<double>("t2b",0));

  mse3RelativePose.get_rotation() = R;
  mse3RelativePose.get_translation() = T;

  double sx = GV3::get<double>("scaleX",1.0);
  double sy = GV3::get<double>("scaleY",1.0);
  double sz = GV3::get<double>("scaleZ",1.0);
  scale[0] = sx; scale[1] = sy; scale[2] = sz;

  if(np==0 || nb==0)
  {
    nbasis=0;
    DBase = NULL;
    return;
  }

  Matrix<> B(nbasis*3,ncomponents);

  //copying the data
  for (unsigned int i=0;i<nbasis*3;i++)
  {
    for(unsigned int j=0;j<ncomponents;j++)
    {
      B[i][j] = 0.0;
    }
  }

  DBase = new Matrix<> (B.T());
}

bool DBases::loadFromFile(string filename)
{
  ifstream ifile;
  ifile.open(filename.c_str(),ifstream::in);

  if(ifile.is_open())
  {
    //preparing temporary data
    vector<vector<double> > tmpdata;
    string tmpline;

    vector<double> tmpvector;
    double tmpdouble;

    while(!ifile.eof() && ifile.good())
    {
      //get the first line from the file. the number of components of each basis is unknown at the begining.
      getline(ifile,tmpline);

      //treat it as a stream from which we are obtaining the first set of data.
      stringstream sstream;
      sstream << tmpline;
      sstream.seekg(0, ios::beg);
      tmpvector.clear();

      //reading each component separatelly
      while(!sstream.eof())
      {
        sstream >> tmpdouble;
        tmpvector.push_back(tmpdouble);
      }

      tmpdata.push_back(tmpvector);
    }

    //flush the last empty vector
    tmpdata.pop_back();

    //now we know the size of the whole matrix
    ncomponents = tmpdata[0].size();
    nbasis = tmpdata.size();

    cout << "Number of basis: " << nbasis << endl;
    cout << "Number of components: " << ncomponents << endl;

    Matrix<> B(nbasis,ncomponents);

    //copying the data
    for (unsigned int i=0;i<nbasis;i++)
    {
      for(unsigned int j=0;j<ncomponents;j++)
      {
        double sc = scale[i%3];
        if(sc!=1.0)
          B[i][j] = tmpdata[i][j]*sc;
        else
          B[i][j] = tmpdata[i][j];
      }
    }

    nbasis /= 3;

    DBase = new Matrix<> (B.T());

    cout << "Deformation base matrix successfuly filled up" << endl;

    ifile.close();

    estimateBaseComponentEnergy();

  }else{
    cerr << "error while trying to load the bases file" << endl;
    return false;
  }
  return true;
}

bool DBases::loadRigidFromFile(string filename)
{
  ifstream ifile;
  ifile.open(filename.c_str(),ifstream::in);

  if(ifile.is_open())
  {
    //preparing temporary data
    vector<vector<double> > tmpdata;
    string tmpline;

    vector<double> tmpvector;
    double tmpdouble;

    int nLines = 0;
    while(!ifile.eof() && ifile.good() && nLines <3)
    {
      //get the first line from the file. the number of components of each basis is unknown at the begining.
      getline(ifile,tmpline);

      //treat it as a stream from which we are obtaining the first set of data.
      stringstream sstream;
      sstream << tmpline;
      sstream.seekg(0, ios::beg);
      tmpvector.clear();

      //reading each component separatelly
      while(!sstream.eof())
      {
        sstream >> tmpdouble;
        tmpvector.push_back(tmpdouble);
      }

      tmpdata.push_back(tmpvector);
      nLines++;
    }

    Matrix<> B(3,ncomponents);

    //copying the data
    for (unsigned int i=0;i<3;i++)
    {
      for(unsigned int j=0;j<ncomponents;j++)
      {
        double sc = scale[i%3];
        if(sc!=1.0)
          B[i][j] = tmpdata[i][j]*sc;
        else
          B[i][j] = tmpdata[i][j];
      }
    }

    mmRigidShape = new Matrix<> (B.T());
    cout << "Rigid Shape matrix successfuly filled up" << endl;

    ifile.close();

    mbRigidShapeLoadedFromFile = true;

  }else{
    cerr << "error while trying to load the bases file" << endl;
    mbRigidShapeLoadedFromFile = false;
    return false;
  }
  return true;
}

//DBases destructor
DBases::~DBases()
{
  delete DBase;
  delete [] mpvBaseComponentEnergy;
  delete [] mpvBaseWeights;
}

void DBases::estimateBaseComponentEnergy(){
  //define the vector of bases energy
  mpvBaseComponentEnergy = new double[nbasis];
  Vector<Dynamic,double,Reference> baseComponentEnergy(mpvBaseComponentEnergy,nbasis);
  baseComponentEnergy = Zeros;
  //Compute each component energy so as to normalize
  Matrix <> B = (*DBase).T();
  for(unsigned int i=0;i<nbasis;i++)
  {
    Vector<> v1 = B[i*3];
    baseComponentEnergy[i] += v1*v1;
    v1 = B[i*3+1];
    baseComponentEnergy[i] += v1*v1;
    v1 = B[i*3+2];
    baseComponentEnergy[i] += v1*v1;
    baseComponentEnergy[i] = sqrt(baseComponentEnergy[i]);
  }
  cout << "baseComponentEnergy " <<  baseComponentEnergy << endl;

  mpvBaseWeights = new double[nbasis];
  Vector<Dynamic,double,Reference> baseWeights(mpvBaseWeights,nbasis);
  double sum=0;
  for(unsigned int i=0;i<nbasis;i++)
  {
    sum += baseComponentEnergy[i];
  }
  baseWeights = baseComponentEnergy;
  sum = 1/sum;
  for(unsigned int i=0;i<nbasis;i++)
  {
    baseWeights[i] *=sum;
  }
  cout << "baseWeights " << baseWeights << endl;
}

void DBases::normalizeBases(bool donorm)
{
  //if not normalized jet
  if(donorm == true && !normalized)
  {
    for(unsigned int i=0;i<nbasis;i++)
    {
      (*DBase).T()[i*3] /= mpvBaseComponentEnergy[i];
      (*DBase).T()[i*3+1] /= mpvBaseComponentEnergy[i];
      (*DBase).T()[i*3+2] /= mpvBaseComponentEnergy[i];
    }
    normalized = true;
    return;
  }

  //if normalized and we want to go back to the standard state
  if(donorm == false && normalized)
  {
    //multiply by the vector of energy pre-calculated
    for(unsigned int i=0;i<nbasis;i++)
    {
      (*DBase).T()[i*3] *= mpvBaseComponentEnergy[i];
      (*DBase).T()[i*3+1] *= mpvBaseComponentEnergy[i];
      (*DBase).T()[i*3+2] *= mpvBaseComponentEnergy[i];
    }
    normalized = false;
    return;
  }
}

void DBases::transformBases()
{
  Matrix<> &b = *DBase;
  for(unsigned int i = 0; i < nbasis; i++)
  {
    for(unsigned int j = 0; j < ncomponents; j++)
    {
        Vector<> v1 = b[j];
        Vector<3> v2 = v1.slice(3*i,3);
        v2 = mse3RelativePose*v2;
        b[j].slice(3*i,3) = v2;
    }
  }
}

void DBases::transformRigid()
{
  Matrix<> &b = *mmRigidShape;
  for(unsigned int i = 0; i < ncomponents; i++)
    b[i].slice<0,3>() = mse3RelativePose*b[i].slice<0,3>();
}
