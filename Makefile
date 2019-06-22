CPP = g++
CPPFLAGSOPENMP = -fopenmp
CPPFLAGSDEBUG = -g -Wall -Wextra
CPPFLAGSSTACK = -Wl,--stack=134217728
# CPPFLAGS = $(CPPFLAGSOPENMP) $(CPPFLAGSSTACK)
CPPFLAGS = $(CPPFLAGSOPENMP) $(CPPFLAGSDEBUG) $(CPPFLAGSSTACK)

default: main.exe

main.exe: main.cpp render.hpp scene.hpp object.hpp ray.hpp vec3.hpp utils.hpp polynome.hpp bezier.hpp Makefile
	$(CPP) $(CPPFLAGS) -o main.exe main.cpp

smallpt.exe: smallpt.cpp Makefile
	$(CPP) $(CPPFLAGS) -o smallpt.exe smallpt.cpp

test.exe: test.cpp Makefile
	$(CPP) $(CPPFLAGS) -o test.exe test.cpp
