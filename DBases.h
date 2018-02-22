// -*- c++ -*-
//Author: Sebastian Bronte <sebastian.bronte@depeca.uah.es>
//Date: 5 May 2011
//DBase.h
//
//DBase class declaration file. Defines principal members and methods to be
//implemented in DBases.cc to load an provide data from the deformation bases.
//


#ifndef _DBASES_
#define _DBASES_

#include <TooN/TooN.h>
#include <TooN/se3.h>
using namespace TooN;

#include <vector>

class DBases
{
  public:
    //constructor and destructor
    DBases(unsigned int np=0, unsigned int nb=0);
    ~DBases();
  
    //specific functions
    //load the whole set of basis from file
    bool loadFromFile(std::string filename);
    //load the rigid component from file, instead of doing it for the map
    bool loadRigidFromFile(std::string filename);

    bool isRigidFromFile(){return mbRigidShapeLoadedFromFile;}
      
    //Get and Set functions
    inline Vector<3> getBasisPoint(const unsigned int bindex,const unsigned int pindex)
    {
        return (DBase->operator[](pindex)).slice(3*bindex,3);
    }

    inline void setBasisPoint(const unsigned int bindex, const unsigned int pindex, Vector<3> &point)
    {
        Matrix <> &B = *DBase;
        Vector<> v1 = B[pindex];
        v1.slice(3*bindex,3) = point;
    }
    //Public Get methods

    //retrieves the number of N Basis in the bases.
    inline unsigned int getNBasis(){return nbasis;}

    //retrieves the components that each basis has got.
    inline unsigned int getNPoints(){return ncomponents;}
    inline double getComponentEnergy(unsigned int index){return (index>nbasis ||index<0) ? 0 : mpvBaseComponentEnergy[index];}

    void normalizeBases(bool donorm);

    SE3<> mse3RelativePose;

    Matrix<> *mmRigidShape;
    
  protected:
    
    //some stuff to make function calls faster
    unsigned int ncomponents;
    unsigned int nbasis;
    
    //Actually this is the main data of this class
    //Matrix associated to each one of the basis of the deformable object
    Matrix<> *DBase;

    double *mpvBaseComponentEnergy;
    double *mpvBaseWeights;
    bool normalized;

    bool mbRigidShapeLoadedFromFile;

    void estimateBaseComponentEnergy();

    void transformBases();
    void transformRigid();

    //scale factors
    Vector<3> scale;
};

#endif
