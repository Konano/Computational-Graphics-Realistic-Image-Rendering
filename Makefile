CPP = g++
CPPFLAGSOPENMP = -fopenmp
CPPFLAGSDEBUG = -g -Wall -Wextra
CPPFLAGSSTACK = -Wl,--stack=134217728
CPPFLAGS = $(CPPFLAGSOPENMP) $(CPPFLAGSDEBUG) $(CPPFLAGSSTACK)

default: main.exe smallpt.exe

main.exe: main.cpp Makefile
	$(CPP) $(CPPFLAGS) -o main.exe main.cpp

smallpt.exe: smallpt.cpp Makefile
	$(CPP) $(CPPFLAGS) -o smallpt.exe smallpt.cpp
