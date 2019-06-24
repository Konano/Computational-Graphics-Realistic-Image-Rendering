CPP = g++
CPPFLAGSOPENMP = -fopenmp
CPPFLAGSDEBUG = -g -Wall -Wextra
CPPFLAGSSTACK = -Wl,--stack=134217728
# CPPFLAGS = $(CPPFLAGSOPENMP) $(CPPFLAGSSTACK)
CPPFLAGS = $(CPPFLAGSOPENMP) $(CPPFLAGSDEBUG) $(CPPFLAGSSTACK)
OPENCV = -I C:\Develop\opencv-3.4.6\build\install\include -L C:\Develop\opencv-3.4.6\build\install\x64\mingw\lib -llibopencv_imgcodecs346 -llibopencv_core346 -llibopencv_highgui346 -llibopencv_imgproc346 -llibopencv_ml346

default: main.exe

main.exe: main.cpp render.hpp scene.hpp object.hpp ray.hpp vec3.hpp utils.hpp polynome.hpp bezier.hpp texture.hpp Makefile
	$(CPP) $(CPPFLAGS) -o mirror.exe main.cpp $(OPENCV)

# main_opencv.exe: main.cpp render.hpp scene.hpp object.hpp ray.hpp vec3.hpp utils.hpp polynome.hpp bezier.hpp Makefile
# 	$(CPP) $(CPPFLAGS) -o main_opencv.exe main.cpp $(OPENCV)

smallpt.exe: smallpt.cpp Makefile
	$(CPP) $(CPPFLAGS) -o smallpt.exe smallpt.cpp

test.exe: test.cpp Makefile
	$(CPP) $(CPPFLAGS) -o test.exe test.cpp
