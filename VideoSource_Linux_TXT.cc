// Author: Sebastián Bronte sebastian.bronte@depeca.uah.es
// Given only the tracks at level 0 from a file, generate a frame with the
// points in. Just display proposals.
#include "VideoSource.h"
#include <cvd/videofilebuffer.h>
#include <cvd/colourspace_convert.h>
#include <cvd/colourspaces.h>
#include <gvars3/instances.h>

using namespace TooN;
using namespace CVD;
using namespace std;
using namespace GVars3;

/**
 * Class constructor
 * @param filename of the file from which the tracks are given
 */
VideoSource_TXT::VideoSource_TXT(string filename)
{
  cout << "  VideoSource_Linux: Opening video source: TEXT..." << endl;

  ImageRef irSize = GV3::get<ImageRef>("VideoSource.Resolution", ImageRef(640,480));

  if(GetTracksFromFile(filename))
      cout << "  ... got text source." << endl;
  else
  {
      cout << "  ... error when loading tracks data" << endl;
      Exceptions::All e;
      e.what = "error trying to load the tracks file";
      throw e;
  }
  mirSize = irSize;
  mpvis = NULL;
  mpcoeffs = NULL;
}

/**
 * Class destructor
 */
VideoSource_TXT::~VideoSource_TXT()
{
  free(mptr);
  free(mpvis);
  free(mpcoeffs);
}

/**
 * Return the size of the image. Not used in this class as it is supposed to be
 * @return size: image size
 */
ImageRef VideoSource_TXT::Size()
{
  return mirSize;
}

/**
 * This function usually retrieves the BW and RGB frames. Not used in this
 * context, since there are no images
 * @param imBW image in black and white. Not available in this case
 * @param imRGB color image in RGB format. Not available in this case
 */
void VideoSource_TXT::GetAndFillFrameBWandRGB(Image<byte> &imBW, Image<Rgb<byte> > &imRGB)
{
  // Nothing to do here ...
}

/**
 * This function returns the tracks from the file name and stores it in the
 * class, as a buffer, to return it later, as they are needed
 * @param filename of the file in which the tracks are stored
 * @return true if the file is correctly loaded, false otherwise
 */
bool VideoSource_TXT::GetTracksFromFile(string filename)
{
  //get the tracks from the file and store the data in the mptr pointer.
  ifstream ifiletracks;
  ifiletracks.open(filename.c_str(),ios::in);
  if(ifiletracks.is_open())
  {
    //preparing temporary data
    vector<vector<double> > tmpdata;
    string tmpline;

    vector<double> tmpvector;
    double tmpdouble;

    while(!ifiletracks.eof() && ifiletracks.good())
    {
      // get the first line from the file. The number of components
      // of each basis is unknown at the begining.
      getline(ifiletracks,tmpline);

      //treat it as a stream from which we obtain the first set of data.
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

    //now we know the size of the whole matrix
    tmpdata.pop_back();
    nPoints = tmpdata[0].size();
    nFrames = tmpdata.size()/2;
    cout << "NTrackPoints: " << nPoints << endl;
    cout << "NFrames: " << nFrames << endl;

    mptr = malloc(sizeof(double)*nFrames*2*nPoints);
    Matrix<Dynamic, Dynamic, double, Reference::RowMajor> Tracks((double*) mptr, nFrames, 2*nPoints);

    //copying the data to the matrix
    for (unsigned int i=0;i<nPoints;i++)
    {
      for(unsigned int j=0;j<nFrames;j++)
      {
        Tracks[j][i*2] = tmpdata[j*2][i];
        Tracks[j][i*2+1] = tmpdata[j*2+1][i];
      }
    }
    ifiletracks.close();
  }
  else
    return false;
  return true;
}

/**
 * This function returns the tracks for the given frame in the parameters
 * @param nframe: frame number from which the tracks are desired
 * @param tracks: the track points for the current frame
 */
void VideoSource_TXT::GetTracksCurrentFrame(unsigned int nframe, vector<Vector <2> > &tracks)
{
  // This matrix is created to encapsulate the data and have an easy access
  // to the data. To speed up, consider to access directly to the elements
  tracks.clear();
  // check if the maximum number of tracks wSas reached
  if(nframe >= nFrames)
  {
    cout << "  ... maximum number of frames on the file reached" << endl;
    Exceptions::All e;
    e.what = "error max number of frames reached";
    throw e;
  }

  Matrix<Dynamic, Dynamic, double, Reference::RowMajor> Tracks((double*)mptr, nFrames, 2*nPoints);

  for(unsigned int i=0;i<nPoints;i++)
  {
    Vector<2> tmpvec;
    tmpvec[0] = Tracks[nframe][i*2];
    tmpvec[1] = Tracks[nframe][i*2+1];
    tracks.push_back(tmpvec);
  }
}

/**
 * Function to retrieve the visibility flag for the current frame
 * @param nframe frame number from which the tracks are desired
 * @param visflags visibility flags for each point in the current frame
 */
void VideoSource_TXT::GetVisibilityCurrentFrame(unsigned int nframe, vector<bool> &visflags)
{
  //visibility is not compulsory to be initialized
  if(mpvis){
    // This matrix is created just to encapsulate the data and have an easy
    // access to the data. To speed up, consider to access directly to the
    // elements
    visflags.clear();
    //check if the maximum number of tracks has reached
    if(nframe >= nFrames)
    {
      cout << "  ... maximum number of frames on the file reached" << endl;
      Exceptions::All e;
      e.what = "error max number of frames reached";
      throw e;
    }

    Matrix<Dynamic, Dynamic, bool, Reference::RowMajor> Visibility((bool*)mpvis, nFrames, nPoints);
    for(unsigned int i=0;i<nPoints;i++)
    {
      visflags.push_back(Visibility[nframe][i]);
    }
  }
}

/**
 * This function loads the visibility flags from the file
 * @param filename which contains the visibility flags
 */
void VideoSource_TXT::GetVisibilityFromFile(string filename)
{
  //allocate memory for mpvis;
  //get the tracks from the file and store the data in the mptr pointer.
  ifstream ifiletracks;
  ifiletracks.open(filename.c_str(),ios::in);
  if(ifiletracks.is_open())
  {
    //preparing temporary data
    vector<vector<bool> > tmpdata;
    string tmpline;

    vector<bool> tmpvector;
    float tmpfloat;
    bool tmpbool;

    while(!ifiletracks.eof() && ifiletracks.good())
    {
      // Get the first line from the file. The number of components of
      // each basis is unknown at the begining.
      getline(ifiletracks,tmpline);

      //treat it as a stream from which we are obtaining the first set of data.
      stringstream sstream;
      sstream << tmpline;
      sstream.seekg(0, ios::beg);
      tmpvector.clear();
      //reading each component separatelly
      while(!sstream.eof())
      {
        sstream >> tmpfloat;
        if(tmpfloat) tmpbool=true;
        else tmpbool=false;
        tmpvector.push_back(tmpbool);
      }

      tmpdata.push_back(tmpvector);
    }

    //now we know the size of the whole matrix
    tmpdata.pop_back();

    //assume that the size of the files are the same as the projection file.
    //further checking will be required

    if(nPoints != tmpdata[0].size() || nFrames < tmpdata.size())
    {
      Exceptions::All e;
      e.what = "error max number of frames reached";
      throw e;
    }

    mpvis = (bool*) malloc(sizeof(bool)*nFrames*nPoints);
    Matrix<Dynamic, Dynamic, bool, Reference::RowMajor> Visibility((bool*) mpvis, nFrames, nPoints);

    //copying the data to the matrix
    for (unsigned int i=0;i<nPoints;i++)
    {
      for(unsigned int j=0;j<nFrames;j++)
      {
        Visibility[j][i] = tmpdata[j][i];
      }
    }
    ifiletracks.close();
  }
  else
  {
    //throw exception
    Exceptions::All e;
    e.what = "error trying to load the track visibility file";
    throw e;
  }
}

/**
 * Loads the Coefficients from the file. It will force to use directly those
 * coefficients instead of computing them
 * @param filename
 */
void VideoSource_TXT::GetCoefficientsFromFile(std::string filename)
{
  //allocate memory for mpcoeffs
  //get the tracks from the file and store the data in the mptr pointer.
  ifstream ifiletracks;
  ifiletracks.open(filename.c_str(),ios::in);
  if(ifiletracks.is_open())
  {
    //preparing temporary data
    vector<vector<double> > tmpdata;
    string tmpline;

    vector<double> tmpvector;
    double tmpdouble;
    //double tmpbool;

    while(!ifiletracks.eof() && ifiletracks.good())
    {
      // Get the first line from the file. the number of components of
      // each basis is unknown at the begining.
      getline(ifiletracks,tmpline);

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


    //now we know the size of the whole matrix
    tmpdata.pop_back();

    nCoeffs = tmpdata[0].size();

    // Assume that the size of the files are the same as the projection file.
    // further checking will be required

    if(nFrames < tmpdata.size())
    {
      Exceptions::All e;
      e.what = "error max number of frames reached";
      throw e;
    }
    //cout << "NTrackPoints: " << nPoints << endl;
    //cout << "NFrames: " << nFrames << endl;

    mpcoeffs = (double*) malloc(sizeof(double)*nFrames*nCoeffs);
    Matrix<Dynamic, Dynamic, double, Reference::RowMajor> Coefficients((double*) mpcoeffs, nFrames, nCoeffs);

    //copying the data to the matrix
    for (unsigned int i=0;i<nCoeffs;i++)
    {
      for(unsigned int j=0;j<nFrames;j++)
      {
        Coefficients[j][i] = tmpdata[j][i];
      }
    }
    ifiletracks.close();
  }
  else
  {
    //throw exception
    Exceptions::All e;
    e.what = "error trying to load the model coefficient file";
    throw e;
  }
}

/**
 * This function retrieves the coefficients for the current frame
 * @param nframe number of frame from which the coefficients are desired
 * @param coeffs the coefficients pre-computed for the current frame
 */

void VideoSource_TXT::GetCoefficientsCurrentFrame(unsigned int nframe, std::vector<double> &coeffs)
{
  // Coefficients are not compulsory to be initialized
  if(mpcoeffs){
    // This matrix is created just to encapsulate the data and have an easy
    // access to the data. To speed up, consider to access directly to the
    // elements
    coeffs.clear();
    //check if the maximum number of tracks has reached
    if(nframe >= nFrames)
    {
      cout << "  ... maximum number of frames on the file reached" << endl;
      Exceptions::All e;
      e.what = "error max number of frames reached";
      throw e;
    }

    Matrix<Dynamic, Dynamic, double, Reference::RowMajor> Coefficients((double*)mpcoeffs, nFrames, nCoeffs);
    for(unsigned int i=0;i<nCoeffs;i++)
    {
      coeffs.push_back(Coefficients[nframe][i]);
    }
  }
}
