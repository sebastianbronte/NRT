# NRT (Real time) Non Rigid Tracking based on PTAM
Source Code Release v1.0

Sebastián Bronte Palacios

Software implementation of the method described in the publications:

 * "Model-Based Real-Time Non-Rigid Tracking" on Sensors, 2017 by S. Bronte,
   L. M. Bergasa, D. Pizarro and R. Barea.
 * "Real-time sequential model-based non-rigid SfM" on the IEEE International
   conference on RObotic Systems (IROS), 2014 by S. Bronte, M. Paladini,
   L. M. Bergasa, L. Agapito and R. Arroyo.
 * "Real-time Sequential non-rigid Structure from Motion" by S. Bronte,
   L. M. Bergasa and D. Pizarro. PhD Thesis. University of Alcalá. Spain,
   July 2017.

This software is based on PTAM (Parallel Tracking and Mapping, as stated
in the paragraphs below) v1.0-r112. To perform the 3D reconstruction, it uses a 
pre-trained model compacted using PCA. From this set of "basis shapes"
and an initial estimate of the pose to start the tracking, after a rigid
camera calibration, the algorithm can start their estimations.

This release is aimed at researchers, developers and other people familiar
with the reconstruction of non-rigid (deformable) surfaces. 

For any request, please e-mail to sebastian dot bronte at gmail.com

As this is based on PTAM software, the following is extracted from its
readme and modified where it is convenient

Parallel Tracking and Mapping for Small AR Workspaces
-----------------------------------------------------

This software is an implementation of the method described in the
paper `Parallel Tracking and Mapping for Small AR Workspaces' by 
Georg Klein and David Murray, which appeared in the proceedings 
of the IEEE/ACM International Symposium on Mixed and Augmented
Reality (ISMAR) 2007.

This release is aimed at experienced software developers and/or
researchers familiar with implementing real-time vision algorithms in
C++ on the platform of their choice.

Questions? E-mail
ptam@robots.ox.ac.uk

REQUIREMENTS
------------
Supported Operating Systems: 
----------------------------
The code was developed on x86/x86-64 Linux. It has also been ported to
Intel-based MacOS X (using X11 for the display). These two operating
systems are supported.

The software can also be compiled on Win32, but this is not our target
platform and so this port is not very clean; see the Windows section later
in this document.

Processor Requirements: 
-----------------------
The software runs at least two processing-intensive threads at the
same time and so needs at least a dual-core machine. Intel Core 2 Duo
processors 2.4GHz+ are fine.

Graphics:
---------
The software requires accelerated OpenGL graphics output. It has been
written and tested only with nVidia cards: this is primarily of
concern under Linux, where the use of an nVidia card and the
proprietary nVidia display drivers are highly recommended. Since the
Linux code compiles directly against the nVidia driver's GL headers,
use of a different GL driver may require some modificºations to the
code.

Video Input:
------------
The software can handle the live input of a video camera with a wide-angle
lens, capable of 640x480x30Hz video capture and an appropriate driver
installation (which is supported by libCVD.) Only monochrome capture is needed
for tracking, colour images are used only for the AR display. A futher
discussion of video input follows later in this file.

Appart from the original input, several new inputs are included, such as:
 * usb camera / webcam
 * 1394 video / firewire
 * recorded video
 * set of text files including tracks, visibility masks, etc.

All the inputs are compiled, so ready to be used through configuration.
In order to select the input, one of the following options on the settings.cfg
file must be included:

 * webcam:
    i wcam
 * firewire:
    i cam
 * recorded video (v4l based):
    i rec
    ifile <videofilename>
 * text:
    i txt
    tracks <tracksfilename>
           (required. Aligned in time 2D point tracks 2FxP,
            F number of frames, P number of points)
    visfile <visibilitymasksfilename>
           (optional. Aligned FxP matrix with a visibility mask:
             1 visible, 0 occluded).
    coeffs <shapecoefficientsfilename>
           (optional. Only for shape debug purposes. The FxK matrix of
            coefficients, K number of meaningfull bases selected, most
            representative ones from PCA. Computed with respect to the current
            set of bases)

Libraries:
----------
The software has three principal dependencies: 

1. TooN - a header library for linear algebra 
2. libCVD - a library for image handling, video capture and computer
vision
3. Gvars3 - a run-time configuration/scripting library, this is a
sub-project of libCVD.

All three above are written by member of the Cambridge Machine
Intelligence lab and are licensed under the LGPL.

Current versions are available from Savannah via CVS:
http://savannah.nongnu.org/projects/toon (for TooN)
http://savannah.nongnu.org/projects/libcvd (for libCVD and GVars)

The latest version of these libraries can be obtained via CVS and ssh:

# export CVS_RSH=ssh
# cvs -z3 -d:pserver:anonymous@cvs.savannah.nongnu.org:/sources/toon co TooN
# cvs -z3 -d:pserver:anonymous@cvs.savannah.nongnu.org:/sources/libcvd co libcvd
# cvs -z3 -d:pserver:anonymous@cvs.savannah.nongnu.org:/sources/libcvd co gvars3

It should be noted, however, that the libraries change rapidly. To
ensure compatibility, it may be useful to download the libraries
corresponding to a time at which they were known to be compatible. To
do so, use the following commands:

# export CVS_RSH=ssh
# cvs -z3 -d:pserver:anonymous@cvs.savannah.nongnu.org:/sources/toon co -D "Mon May 11 16:29:26 BST 2009" TooN
# cvs -z3 -d:pserver:anonymous@cvs.savannah.nongnu.org:/sources/libcvd co -D "Mon May 11 16:29:26 BST 2009" libcvd
# cvs -z3 -d:pserver:anonymous@cvs.savannah.nongnu.org:/sources/libcvd co -D "Mon May 11 16:29:26 BST 2009" gvars3

The installation of these libraries is described below.

Additionally, the following libraried must be installed:
4. CGAL 3.9
5. OpenMP
6. OpenCV 3.1.1 or greater, including the 3rd party package that allows some additional feature algorithms such as SURF and SIFT.

INSTALLATION
------------
Installation of the dependencies
--------------------------------

The three dependent libraries are all compiled and installed using the
familiar ./configure; make; make install system, however the following
points are worth noting.

On Linux, the following libraries (and their -devel versions) are
required: blas, lapack, perhaps libgfortran, ncurses and libreadline
(optional, for GVars3), libdc1394 (and maybe libraw1394)
for firewire capture, optionally libtiff, libjpeg, libpng.

(On OSX, these or their equivalents should all be available with the
system.)

The order of installation should be 
1. TooN, 2. libCVD, 3. GVars3;

TooN installation is trivial, since it's only a bunch of headers.

For libCVD, I recommend the following configure options:
# export CXXFLAGS=-D_REENTRANT
# ./configure --without-ffmpeg

In case the video input (V4L2) is required, the following configuration flags
must be provided:
# ./configure --enable-gpl --with-ffmpeg --with-v4l2buffer

if the compiler hangs on one of the FAST detectors these can be
disabled with configure --disable-fast7 (for example.)
Documentation can be generated (if doxygen is installed) by running
# make docs

For GVars3, I recommend the following configure options:
# ./configure --disable-widgets

After that, download the compressed file of CGAL3.9 from the internet, unzip
and follow the instructions inside the package to perform the installation.

Finally, install OpenCV 3.1.1 or greater, including the 3rd party package
to try with the OpenCV implementation of SIFT and SURF detectors. To do that,
follow the instructions on the package and perform an specific installation
using ccmake . and cmake ., to properly configure and install the package.

Compiling the Software
----------------------
The source code is within the main directory.

The Makefile can then be edited to reference any custom include or linker
paths which might be necessary (depending on where the dependencies were
installed.)

By default, several video sources are compiled and chosen by configuration,
instead of compilation time. Video recording option is also available, with
the proper configuration.

The second step, for Linux, is to set up the correct video source. Two
files are provided, VideoSource_Linux_DV.cc and
VideoSource_Linux_V4L.cc, which work with the Unibrain Fire-i and the
Logitech Quickcam Pro 5000 respectively. The DV version is compiled by
default; edit the Makefile to switch to the V4L version instead, if
needed. Other cameras may require manual editing of the video input
files, e.g. to change the videobuffer's colourspace.

Other video source classes are available with libCVD. Finally, if a
custom video source not supported by libCVD is required, the code for
it will have to be put into some VideoSource_XYZ.cc file (the
interface for this file is very simple.)

The software can then be compiled with the command 
# make

This builds two target executables: NRT and CameraCalibrator.

RUNNING THE SOFTWARE
--------------------
Calibrating the Camera
----------------------
CameraCalibrator should be run first to obtain a camera calibration
(and to verify that video input is in fact working.) This requires the
user to point the camera at a checker-board calibration pattern; any
checkerboard of any size will do, a sample is included as
calib_pattern.pdf.

The camera calibrator attempts to find square corners in the image,
and then to link these together. This is indicated by fragments of
squares appearing in the image. The calibrator is not very robust; If
no squares are detected, some things to try would be:
- Modify camera settings to remove any sharpening; Sharpening
  artefacts break the calibrator, and also reduce performance of the
  tracker.
- Adjust the calibrator's settings, for example increase the value of
  CameraCalibrator.BlurSigma

When the camera is in a pose in which some portions of the grid are
detected, the user should press the `GrabFrame' button to add a
snapshot of that frame to the optimiser. After a few frames from
different poses have been added, pressing `Optimize' button
iteratively calculates the camera parameters. When the user is
satisfied with convergence (the RMS error should be no more than
around 0.3 pixels) pressing `Save' stores the camera calibration in a
file camera.cfg.

For sinthetically obtained dataset, this step could be skipped by directly
creating a custom camera.cfg file including the normalized (w.r.t the size of
the video input) focus in u and v coordinates, optical center in u and v and
the normalized radial distortion coefficient.

Running the Tracker
-------------------
Once a calibration has been stored, invoking NRT runs the
tracker.

This step requires a setup as the proper initialization has not been completed.
The inputs needed to start the sequence are the following and should be
included when the program starts:

 * r00,r01,r02,r10,r11,r12,r20,r21,r22:
     The initial rotation matrix approximation (row,col).
 * t0,t1,t2:
     The initial translation vector approximation.
 * bases <basesPointWiseFile>:
     The deformation bases as a set of basis shapes, result from the
     computation as depicted in the related publications, shaped in 3KxP
 * rigid <rigidPointWiseFile>
     The estimation of the rigid shape of the object to be tracked.
 * FeatureType <featuretype>
     This type could be one of the following: SIFT, SURF, PTAM, ORB, BRISK, KAZE,
     AKAZE. PTAM type requieres the compilation and substitution of certain parts
     of the code, as it is not compatible with the OpenCV parts, but it will be
     conveniently integrated.
 * MatchingType <matchingtype>
     This type could be one of the following: BruteForce, BruteForce_L1,
     BruteForce_L2, BruteForce_Hamming
 * roiu1, roiv1, ... roiu4, roiv4:
     The ROI (pixel coordinates u,v) in which the target object is located.

The former PTAM initialization is not available, as it is thought for
near-planar rigid objects. This algorithm uses the input data and a region of
interest (if needed) in which the object is located on the first frame.

General Use of GVars
--------------------
Both programs rely on GVars for a console user-interface. In the
terminal window, a user may inspect a list of all tweakable variables
by typing 
> gvarlist 
(or gvarlist -a for a more complete list) and then modify variables by
typing
> Variable_Name = new_value
This allows the user to tweak a number of the program's features; for
example, if the quality of video input is dubious, typing
DrawFASTCorners=1 in the tracker allows the user to inspect the
response of the FAST corner detector.

VIDEO SOURCES
-------------

This version of deformable PTAM based tracker is mainly tested with text-based
and pre-recorded sequences. Although the live video input options are
available, the algorithm is not tested with them, as it is hard to give a
priory estimation of the rigid shape that matches the pre-loaded deformation
model, as well as the initial pose.

The rigid counterpart of the software was developed with a Unibrain Fire-i
colour camera, using a 2.1mm M12 (board-mount) wide-angle lens. It also runs
well with a Logitech Quickcam Pro 5000 camera, modified to use the same 2.1mm
M12 lens.

Wide-angle lenses are important for natural feature trackers and SLAM
systems, and the software's performance using a zoomier lens is likely
to be compromised. It will require more key-frames, robustness to
rotation will decrease, and relocaliser performance will drop.

How wide is wide-angle? The first number in camera.cfg is a normalized
horizontal focal length (fx): for our wide-angle lenses fx=0.57. Up to
fx<=1.0 the system would be expected to run fine, beyond that
performance will likely be degraded.

Independent of OS, it is important to turn down in-camera sharpening to the
point that no sharpening artifacts appear! Sharpening artifacts produce
spurious corners from image noise, move the location of real corners, and
break scale and rotation invariance. Expect poor performance if running on
an over-sharpened image. For example, on the unibrain fire-i, turn sharpness
down from the default 80 to 25.

Linux fire-i notes: Firewire with libCVD has been tested against the
old-style (pre-juju) firewire stack, using either the libDC-1 and
libDC-2 series. If your distribution ships with firewire support in
the form of the experimental juju stack (e.g. Fedora 8 and 9) you may
experience video input lockups, you should install packages for the
old firewire stack instead.

Linux Logitech Quickcam notes: the Logitech Quickcam pro 5000 is
supported by the linux-uvc driver and can be used with a
CVD::V4LBuffer<yuv422>.

MacOS X notes: Video input properties are reset every time the software
is run, and so settings will have to be adjusted every time
(e.g. turning down in-camera sharpening.) 

Windows notes: If attempting to port to windows, note that we have
experienced poor performance in conjunction with DSVL: this seems to be a
threading issue, as the video source also uses a thread, which can get
starved.  Using the unibrain or point grey APIs as appropriate would be
preferable. Note also that YUV->RGB->greyscale produces notable artefacts as
opposed to direct YUV->greyscale. At the moment, only a CMU1394 interface is
included, this works fine.

There are a couple of input sources extra provided to be executed for the
tests, which are pre-recorded video input and text-based input.

The pre-recorded video processing is given by V4Linux drivers, in which the
file VideoSource_Linux_V4L.cc/h is implementing an interface that provides
frames to the system.

Regarding the text input, several files are expected, as already indicated
in the video input part of the requirements.

RESULTS COMPILATION: VIDEO RECORDING AND RESULT TEXT FILES GENERATION
---------------------------------------------------------------------

The output video of the processing can be recorded if, in the settings file the following flags are indicated (as previously provided with the code within PTAMM, PTAM for multiple mapping):

 * SaveFIFO=1

The compiled output files are automatically generated following the pattern given by the following configuration option:

 * odata=<stringPattern>

From the depicted string pattern, the following files will be generated:
 * <stringPattern>_features.txt
   This file is structured as follows. First a row containing the number of
   expected features matched with the provided model for the current frame.
   Then, the following rows contain the 2D coordinates where the feature is
   detected, its index, its computed error w.r.t the model, its visibility
   flag, its level in the feature pyramid, whether it is an outlier, and the
   3D computed coordinates for that point (w.r.t world coordinates).

 * <stringPattern>_tracking_poses.txt
   This file contains the computed pose of the camera w.r.t world coordinates
   for each frame. First the rotation coordinates (r0, r1, r2) and then the
   translation coordinates (t0, t1, t2).

 * <stringPattern>_quality.txt
   The detected tracking quality for each frame.

 * <stringPattern>_rawtracks.txt
   Similarly to the _features.txt file, this file stores the whole set of
   detected features on the current frame. In the _features file the only
   features saved are the ones matched with the model. In this file the
   total amount of 2D features (one per row) are stored.

 * <stringPattern>_coefs.txt
   In this file, the shape deformation coefficients are stored, one set of
   coefficients per row.

Activating these data output functions, as one could expect, lower the
performance of the overall algorithm.

DEFORMABLE MODEL LOADING
------------------------

To start the processing of the sequence, whether it was pre-recorded, live,
or text-based, a point-wise (sparse) model must be included to start with
the estimations.

This model should be computed as indicated in the papers, but as a summary,
a previously given set of points over time can be provided as a starting point.
Then PCA is computed on that set and the most significant 

THE SOURCE CODE 
--------------- 
The main documentation of the code is in the code's comments. What
follows here is an overview of the tracking system.
1
The tracker runs from main.cc, which spawns a System-class which
starts two main threads:

1. The Tracker-class thread, driven by the System class event loop,
which performs real-time tracking and video I/O and graphical display;

2. The MapMaker-class thread, which updates the map in the
background. On this deformable version of the tracking, this thread
is disabled, as the model loading overrides this capacity.

Both threads access a common Map-class data structure. The Map-class
contains two main arrays: a vector of MapPoint-structs, and a vector
of KeyFrame-structs. Together these make up the map.

For the case of deformable objects, this Map-class data structure has been
extended, so as to incorporate the multiple shapes of the linear modelling
used, the coefficients, etc.

The bulk of the functionality is contained within the classes Tracker,
NRTracker and VisualNRTracker, which also make use of auxiliary classes. The
tracker directly inherited from PTAM uses the Relocaliser-class to recover
from failure, for more advanced versions of feature-based matching approaches,
This class is recommended not to be used. Many bits of code use
ATANCamera-class, which is a model of the camera's intrinsic parameters.

The MapMaker class has been removed as in this version of the project
tracking is the main part.

COMPILING AND RUNNING ON WINDOWS (from PTAM)
--------------------------------
The software compiles fine on Windows, but the software remains a console
application (gvars interface is in a terminal window) and we experience some
frame-drop issues which don't seem to arise on other platforms.

To compile (and maybe also to run) the software on windows, the following
libraries are needed:

Lapack and BLAS - available pre-compiled e.g. from 
http://www.fi.muni.cz/~xsvobod2/misc/lapack/

pthreads from 
http://sourceware.org/pthreads-win32/

GLEW from 
http://glew.sourceforge.net/

CMU1394 camera driver from 
http://www.cs.cmu.edu/~iwan/1394/

libjpeg for win32 e.g. from 
http://gnuwin32.sourceforge.net/packages/jpeg.htm

Obtain LibCVD, TooN, and GVars3 from CVS as above. To compile and install,
libCVD and GVars3 come with msdev project files, TooN is just a bunch of
headers (copy it to your include directory.)

Copy the contents of PTAM/Build/Win32 to PTAM; then open
PTAM.sln, which contains projects to build the Camera Calibrator and
the Tracker. Edit include and lib directories as appropriate. Good luck!

The previous steps come from the original version of PTAM. As there are new
libraries, such as CGAL3.9 and OpenCV, required for this to work, they should
be properly installed as appropriate in Windows.
