# Detect system OS.
ifeq ($(OS),Windows_NT)
    detected_OS := Windows
else
    detected_OS := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif

CPP = g++
CPPFLAGSOPENMP = -fopenmp
CPPFLAGSC11    = -std=c++11
CPPFLAGSGDB    = -g
CPPFLAGSO3     = -O3
CPPFLAGSDEBUG  = -Wall -Wextra
CPPFLAGSSTACK  = -Wl,--stack=134217728
OPENCV = -I C:\Develop\opencv-3.4.6\build\install\include -L C:\Develop\opencv-3.4.6\build\install\x64\mingw\lib -llibopencv_imgcodecs346 -llibopencv_core346 -llibopencv_highgui346 -llibopencv_imgproc346 -llibopencv_ml346

ifeq ($(detected_OS),Windows)
	CPPFLAGS = $(CPPFLAGSOPENMP) $(CPPFLAGSGDB) $(CPPFLAGSDEBUG) $(CPPFLAGSSTACK) $(OPENCV)
	TARGET = .exe
else
	CPPFLAGS = $(CPPFLAGSOPENMP) $(CPPFLAGSGDB) $(CPPFLAGSC11)
	TARGET =
endif

default: main.exe

test: main.cpp render.hpp scene.hpp object.hpp ray.hpp vec3.hpp utils.hpp polynome.hpp bezier.hpp texture.hpp hitpoint.hpp Makefile
	$(CPP) -o main$(TARGET) main.cpp $(CPPFLAGS)
	gdb main

main.exe: main.cpp render.hpp scene.hpp object.hpp ray.hpp vec3.hpp utils.hpp polynome.hpp bezier.hpp texture.hpp hitpoint.hpp Makefile
	$(CPP) -o main$(TARGET) main.cpp $(CPPFLAGS)

