/*
      C++ program to study fluids.
      by reading in instantaneous configurations 
      and creating neighbor gx(r)
      */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)<(B)?(A):(B))
#define anint(A) ((int)((A)+((A)>=0?.5:-.5)))/* stupid c doesn't know anint */

const int DIM=3;
const int MAXLINELEN=256;
const int NSAVE=10;      /* Save files every NSAVE configs read in */

#include "flags.H"
#include "histogram.h++"
#include "gor_io.h++"

/*----------------------------------------------------------*/
void gofr(int natoms,double x[],double cell,HISTOGRAM &hist)
{
  int i,j,k;
  double dx[DIM],r,rcell;

  rcell = 1./cell;
  double gmax = hist.xmax();
  double gmin = hist.xmin();

  for (i=0;i<natoms;i++){
    for (j=0;j<natoms;j++){
      if(i != j ){ /* don't do same atom */
	for(k=0;k<DIM;k++) dx[k]=x[i*DIM+k]-x[j*DIM+k];
	for(k=0;k<DIM;k++) dx[k] -= anint(dx[k]*rcell)*cell;
	r = 0.;
	for(k=0;k<DIM;k++) r += dx[k]*dx[k];
	r = sqrt(r);
	if (r<gmax && r>gmin) hist.bin(r);
      }
    }
  }
}

/*----------------------------------------------------------*/
void save(int natoms,int nsave,double avvol,HISTOGRAM &hist)
{
  int j;
  FILE *fp;

  fprintf(stdout,"Saving files (%s) with %d configs read in.\n",
	  hist.filename(),nsave);

  double dx = hist.dx();
  double fr=avvol/(4.*M_PI*dx*(double)(natoms*natoms*nsave));

  if((fp = fopen(hist.filename(),"w"))==NULL){
    fprintf(stderr,"Error: Can't open %s\n",hist.filename());
    exit(1);
  }
  fprintf(fp,"# %d %d %g\n",natoms,nsave,cbrt(avvol));
  double *grt = hist.hist();
  double gmin = hist.xmin();
  int npt = hist.ntp();
  for(j=0;j<npt-1;j++){
    double pt = (j+.5)*dx+gmin;
    fprintf(fp,"%lg %lg\n",pt,grt[j]*fr/(pt*pt));
  }
  fclose(fp);
}

/* ------------------------------------------------------------------*/
int main(int argc,char **argv)  /* pass in command line arguments */
{
  int nsave,natoms,nconf;
  double *x,dt,cell,avvol;
  FLAGS flags;
  HISTOGRAM hist;
  FILE *fp;

  /* read in command line arguments */
  flags.read(argc,argv);
  hist.set(flags.n(),flags.xmin(),flags.xmax(),flags.outfile());
  
  fp = readHeader(flags.infile(),natoms,nconf,dt);

  if(!flags.q()){
    fprintf(stdout,"# of atoms = %d, ",natoms);
    fprintf(stdout,"steps per # of configurations = %d\n",nconf);
    fprintf(stdout,"dt = %lg ps/time step\n",dt);
  }
  /* allocate memory */

  x = (double *)malloc(DIM*natoms*sizeof(double));
  avvol = 0;
  nsave = 0;
  while(getConfig(fp,natoms,x,cell)){
    avvol += cell*cell*cell;
    ++nsave;
    gofr(natoms,x,cell,hist);

    if(!(nsave%NSAVE)){save(natoms,nsave,avvol/nsave,hist);}
  }
  fclose(fp);

  if((nsave%NSAVE)){save(natoms,nsave,avvol/nsave,hist);}

  fprintf(stdout,"Done with %d configs averaged\n",nsave);

  return 0;
}

