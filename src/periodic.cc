/*! \file periodic.cc
  file the contains the function that pertain to the periodic nature
  or boundry conditions of the simulation.
*/

#include "periodic.h++"

// Static variables
static int iperd=-1; //!< the periodicity of the simulation
static double box[DIM]; //!< the box length
static double rbox[DIM]; //!< the reciprocal of the box length
static double rad2; //!< the square of the radius used in spherical boundry conditions

//! set the periodicity and boxlength in periodic
void setIbox(int iperdt,double boxt){
  int i;
  iperd = iperdt;
  for(i=0;i<DIM;++i){ box[i] = boxt;}
  for(i=0;i<DIM;++i){ rbox[i] = 1./box[i];}
  rad2 = box[0]*box[0];
}

int getIperd(void){return iperd;} //!< returns the periodicity
double getBox(void){return box[0];} //!< returns the cell length
double getVol(void){return box[0]*box[1]*box[2];} //!< return the volume

/*! 
  scale coords to so they are fractions of the cell which
  will now have a length of 1
*/
void scalecoord(int nintr,double xdis[]){
  int i,j;
  nintr *= DIM; 
  for(i=0;i<nintr;i+=DIM){
    for(j=0;j<DIM;++j){
      xdis[i+j] *= rbox[j];
    }
  }
}

/*! unscale coords to original values */
void unscalecoord(int nintr,double xdis[]){
  int i,j;nintr *= DIM; 
  for(i=0;i<nintr;i+=DIM){
    for(j=0;j<DIM;++j){
      xdis[i+j] *= box[j];
    }
  }
}

//! assume xdis already scaled and just put back into cell
void addcell(int nintr,double xdis[],int icell[]){
  int i,j;
  nintr *= DIM; 
  for(i=0;i<nintr;i+=DIM){
    for(j=0;j<DIM;++j){
      icell[i+j] += anint(xdis[j+i]);
      xdis[j+i] -= anint(xdis[j+i]);
    }
  }
}

void scaletocell(int nintr,double xdis[],int icell[]){
  scalecoord(nintr,xdis);
  addcell(nintr,xdis,icell);
}

void scaletoreal(int nintr,double xdis[],int icell[]){
  int i,j;
  nintr *= DIM; 
  for(i=0;i<nintr;i+=DIM){
    for(j=0;j<DIM;++j){
      xdis[j+i] = (xdis[j+i]+icell[j+i])*box[j]; 
    }
  }
}

/*! 
  applies the spherical boundary conditions
  \param nintr is the number distances to scale
  \param xdis is a double array of distances which will be returned
  \param vdis is a double array of velocities to reflect 
  with the boundary conditions applied.
*/
void spherical(int nintr,double xdis[],double vdis[])
{
  int i,j;
  if(iperd != 4){
    // ERROR();  // -2 r . v / |r| + v
  }
  for(i=0;i<nintr;++i) {
    double r2 = 0;
    int ii = i*DIM;
    for(j=0;j<DIM;++j){r2 += xdis[ii+j]*xdis[ii+j];}
    if(r2 > rad2){
      // put particle at boundry
      double fact = sqrt(rad2/r2);
      for(j=0;j<DIM;++j){ xdis[ii+j] *= fact;}
      // bounce velocity 
      double proj = 0.0;
      for(j=0;j<DIM;++j) proj += vdis[ii+j]*xdis[ii+j];
      //#ifdef DEBUG
      if(proj < 0 ) {printf("Wrong bounce in sherical periodic conditions!\n");}
      //#endif
      for(j=0;j<DIM;++j) vdis[ii+j] += -2.*proj*xdis[ii+j]/rad2;
    }
  }
}

/*! 
  applies the periodic images to the vector xdis
  \param nintr is the number distances to scale
  \param xdis is a double array of distances which will be returned
  with the periodic boundry conditions applied.
*/
void periodic(int nintr,double xdis[])
{
  int i;
  nintr *= DIM;

  switch (iperd){ 
  case 4: // spherical boundry conditions don't move particles, all done in spherical after each step)
    break;
    // cartesian periodic boundary conditions fall through
  case 3: for(i=0;i<nintr;i+=DIM) xdis[2+i] -= box[2]*anint(xdis[2+i]*rbox[2]);
  case 2: for(i=0;i<nintr;i+=DIM) xdis[1+i] -= box[1]*anint(xdis[1+i]*rbox[1]);
  case 1: for(i=0;i<nintr;i+=DIM) xdis[i]   -= box[0]*anint(xdis[  i]*rbox[0]);
  case 0: break;
  default: fprintf(stderr,"%d periodic not allowed!\n",iperd);exit(1);break;
  }
#ifdef DEBUG
  for(i=0;i<nintr;i+=DIM) {
    double r2 = 0;
    for(int j=0;j<DIM;++j){r2 += xdis[i+j]*xdis[i+j];}
    printf("periodic: i=%d x=%g y=%g z=%g r=%g\n",i,xdis[i],xdis[i+1],xdis[i+2],sqrt(r2));
  }
#endif
}

/*! 
  applies the periodic images to the vector xdis assuming the values
  are already scaled so that the box lenght is one.
  \param nintr is the number distances to scale
  \param xdis is a double array of distances which will be returned
  with the periodic boundry conditions applied.
*/
void iperiodic(int nintr,double xdis[])
{
  int i,j;
  nintr *= DIM;

  for(i=0;i<nintr;i+=DIM) for(j=0;j<DIM;++j) xdis[i+j] *= rbox[j];

  switch (iperd){ // fall through
  case 3: for(i=0;i<nintr;i+=DIM) xdis[2+i] -= anint(xdis[2+i]);
  case 2: for(i=0;i<nintr;i+=DIM) xdis[1+i] -= anint(xdis[1+i]);
  case 1: for(i=0;i<nintr;i+=DIM) xdis[  i] -= anint(xdis[  i]);
  case 0: break;
  default: fprintf(stderr,"%d periodic not allowed!\n",iperd);exit(1);break;
  }
}

void resetcoords(int nintr,double xdis[],double vdis[]){
  if(iperd==4){
    spherical(nintr,xdis,vdis);
  }
}

