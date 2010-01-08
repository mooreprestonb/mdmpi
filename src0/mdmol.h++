/*! \file mdmol.h++
  file to hold the prototypes and meta includes for 
  the program mdmol
*/

#ifndef _MDMOL_H_
#define _MDMOL_H_

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#ifdef PARA
#include <mpi.h>
#endif

#include "flags.H"
#include "dmatrix.h++"
#include "simvars.h++"
#include "coords.h++"
#include "energy.h++"
#include "neighbor.h++"
#include "periodic.h++"
#include "parallel.h++"
#include "timing.h++"
#include "integrate.h++"
#include "lj_energy.h++"
#include "electrostatic.h++"
#include "freeze.h++"
#include "mdOut.h++"
#include "mathutil.h++"
#include "zerotm.h++"
#include "mdinclude.h++"
void *cmalloc(size_t);

void setup(FLAGS &,SIMVARS &,COORDS &,NEIGHBOR &,ENERGY &);

void freeEnergyScr(void);
double pot(COORDS &coords,NEIGHBOR &ngbr,int);
void force(COORDS &coords,NEIGHBOR &ngbr);
double getpot(double *,double *,double *,double *,double *,double *);
double getke(NEIGHBOR &);

void getMyvals(int npart,double x[],double xp[],int ipmap[]);
void scatMyvals(int npart,double x[],double xp[],int ipmap[]);
int noReq(int rank);
void doneReq(int nproc[]);
int packxi(int,int,int[],int[],double xib[] ,double xjb[],int[],double[]);
int alldone(int nr,int iproc[],int n);

/*---------------------------------------------------------------------*/

//#define vlrc0 0. 
//#define wlrc0 0.

enum TAGS {TAGREQ,TAGSERVI,TAGSERVD,TAGUPDI,TAGUPDD,TAGCMI,TAGCMD,DONETAG};
#define ITERMAX 100000

/*---------------------------------------------------------------------*/

#endif
