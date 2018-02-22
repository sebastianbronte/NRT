// Author: Sebastián Bronte sebastian.bronte@depeca.uah.es
// In order to read data from a saved video sequence
#include "VideoSource.h"
#include <cvd/videofilebuffer.h>
#include <cvd/colourspace_convert.h>
#include <cvd/colourspaces.h>
#include <gvars3/instances.h>

using namespace TooN;
using namespace CVD;
using namespace std;
using namespace GVars3;

VideoSource_FF::VideoSource_FF(string FileName)
{
  cout << "  VideoSource_Linux: Opening video source: FFMPEG..." << endl;
  VideoFileBuffer<Rgb<byte> >* pvb = new VideoFileBuffer<Rgb<byte> >(FileName);
  pvb->on_end_of_buffer(CVD::VideoBufferFlags::UnsetPending);
  mirSize = pvb->size();
  mptr = pvb;
  cout << "  ... got video source." << endl;
}

VideoSource_FF::~VideoSource_FF()
{
  delete (VideoFileBuffer<Rgb<byte> >*) mptr;
}

ImageRef VideoSource_FF::Size()
{
  return mirSize;
}

void VideoSource_FF::GetAndFillFrameBWandRGB(Image<byte> &imBW, Image<Rgb<byte> > &imRGB)
{
  VideoFileBuffer<Rgb<byte> >* pvb = (VideoFileBuffer<Rgb<byte> >*) mptr;
  VideoFrame<Rgb<byte> > *pVidFrame = pvb->get_frame();
  convert_image(*pVidFrame, imBW);
  imRGB.copy_from(*pVidFrame);
  pvb->put_frame(pVidFrame);
}

/**
 * Loads the Coefficients from the file. It will force to use directly those
 * coefficients instead of computing them
 * @param filename
 */
void VideoSource_FF::GetCoefficientsFromFile(std::string filename)
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

    while(!ifiletracks.eof() && ifiletracks.good())
    {
      //get the first line from the file. the number of components of each basis is unknown at the begining.
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
    nFrames = tmpdata.size();
    nCoeffs = tmpdata[0].size();

    //assume that the size of the files are the same as the projection file.
    //further checking will be required
    if(nFrames < tmpdata.size())
    {
      Exceptions::All e;
      e.what = "error max number of frames reached";
      throw e;
    }

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

void VideoSource_FF::GetCoefficientsCurrentFrame(unsigned int nframe, std::vector<double> &coeffs)
{
  //coefficients are not compulsory to be initialized
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
