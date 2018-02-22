# DO NOT DELETE THIS LINE -- make depend depends on it.
#MY_CUSTOM_INCLUDE_PATH = ../utils
#MY_CUSTOM_LINK_PATH = ../utils


# Edit the lines below to point to any needed include and link paths
# Or to change the compiler's optimization flags
MY_CUSTOM_INCLUDE_PATH=/homes/sbronte/PTAM_old/include
MY_CUSTOM_LINK_PATH=
CC=g++

COMPILEFLAGS=-I $(MY_CUSTOM_INCLUDE_PATH) -D_LINUX -D_REENTRANT -Wall -march=nocona -msse3 -lCGAL -frounding-math -g #-O3 -fopenmp
LINKFLAGS=-L $(MY_CUSTOM_LINK_PATH) -lmpfr -lgmp -lGVars3 -lcvd -llapack -ldc1394 -lraw1394 -lpng -ljpeg -lavcodec -lavformat -lswscale -ltiff -lGL -lCGAL -frounding-math -lblas -lopencv_core -lopencv_highgui -lopencv_calib3d -lopencv_imgproc -lopencv_features2d -lopencv_xfeatures2d -g #-O3 -lgomp -fopenmp

VIDEOSOURCE_LIVE = VideoSource_Linux_DV.o #input from firewire camera
VIDEOSOURCE_REC = VideoSource_Linux_FF.o #input from recorded video sequence. ffmpeg needed to be installed correctly
VIDEOSOURCE_WCAM = VideoSource_Linux_V4L.o #input from webcam.
VIDEOSOURCE_TXT = VideoSource_Linux_TXT.o #input from raw text files containing the tracks, visibility masks, etc.

OBJECTS=	main.o\
		GLWindow2.o\
		GLWindowMenu.o\
		$(VIDEOSOURCE_LIVE)\
		$(VIDEOSOURCE_REC)\
                $(VIDEOSOURCE_WCAM)\
		$(VIDEOSOURCE_TXT)\
		System.o\
		ATANCamera.o\
		KeyFrame.o\
		MapPoint.o\
		Map.o\
		SmallBlurryImage.o\
		ShiTomasi.o \
		HomographyInit.o \
		PatchFinder.o\
		Relocaliser.o\
		MiniPatch.o\
		Tracker.o\
		MapViewer.o\
		TrackWriter.o\
		DBases.o\
		NRTracker.o\
	        VisualNRTracker.o\
		VisualNRTrackerV1.o\
	        EstimationUtils.o\
	        Interpolation.o\
		Interpolation_LN.o\
	        DeformablePatch.o\
                OpenCVUtilities.o\
		DescriptorMatcherOCV.o

CALIB_OBJECTS=	GLWindow2.o\
		GLWindowMenu.o\
		$(VIDEOSOURCE_LIVE)\
		$(VIDEOSOURCE_REC)\
		$(VIDEOSOURCE_WCAM)\
		CalibImage.o \
		CalibCornerPatch.o\
		ATANCamera.o \
		CameraCalibrator.o

All: NRT CameraCalibrator

NRT: $(OBJECTS)
	$(CC) -o NRT $(OBJECTS) $(LINKFLAGS)

CameraCalibrator:$(CALIB_OBJECTS)
	$(CC) -o CameraCalibrator $(CALIB_OBJECTS) $(LINKFLAGS)

%.o: %.cc
	$(CC) $< -o $@ -c $(COMPILEFLAGS)

clean:
	rm *.o

depend:
	rm dependecies; touch dependencies
	makedepend -fdependencies $(INCLUDEFLAGS) $(MOREINCS) *.cc *.h

-include dependencies









