
SHELL = /bin/sh

CXX = g++
MV = /bin/mv
CP = /bin/cp
RM = /bin/rm 
TOUCH = /usr/bin/X11/touch
LIBS =  -lm
DEST = $(HOME)/bin/
CXXFLAGS = -g -O2 -Wall


mdmol :
	$(MAKE) -C src

gor :
	$(MAKE) -C analysis

clean :
	$(MAKE) -C src clean
	$(MAKE) -C analysis clean
	$(RM) -rf mdmol
