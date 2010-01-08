#  Makefile

# common definitions...
# HOME = /home/moore

# Compiler flages, paths for include files and libraries
# CC = insure -n32
# CFLAGS = -g -fullwarn -DPARA
# LIBS = -lmpi -lm 
# CC = cc -n32 -fullwarn -g
# CFLAGS = -O -fullwarn
# INCLUDE = -I/home/moore/software/fftw-1.3/src
# LLIBS = -L/home/moore/software/fftw-1.3/src

CXX = g++
MV = /bin/mv
CP = /bin/cp
RM = /bin/rm 
TOUCH = /usr/bin/X11/touch
LIBS =  -lm
DEST = $(HOME)/bin/
CXXFLAGS = -g -O2 -Wall

#targets

all :: gor gzero gz

gz : gz.cc histogram.o gor_io.o
	$(CXX) -o $@ $(CXXFLAGS) $(INCLUDE) gz.cc histogram.o gor_io.o $(LLIBS) $(LIBS)

gor : gor.cc histogram.o gor_io.o
	$(CXX) -o $@ $(CXXFLAGS) $(INCLUDE) gor.cc histogram.o gor_io.o $(LLIBS) $(LIBS)

gzero : gzero.cc histogram.o gor_io.o
	$(CXX) -o $@ $(CXXFLAGS) $(INCLUDE) gzero.cc histogram.o gor_io.o $(LLIBS) $(LIBS)

histogram.o : histogram.cc histogram.h++
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) histogram.cc $(LLIBS) $(LIBS)

gor_io.o : gor_io.cc gor_io.h++
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) gor_io.cc $(LLIBS) $(LIBS)

clean:
	$(RM) -rf gor_io.o histogram.o gz gor gzero
