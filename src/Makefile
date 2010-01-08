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

all :: mdmol

include Makefile.mdmol

mdmol : $(MDFILES) mdmol.h++ 
	$(CXX) -o $@ $(CXXFLAGS) $(INCLUDE) $(MDFILES) $(LLIBS) $(LIBS)
	$(MV) $@ ../$@

