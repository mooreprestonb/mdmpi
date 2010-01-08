/*! \file force.cc 
  file to calculate the forces.
  this just collects the particles that interact
  via the neighbor lists and then call 
  functions to actually calculate the force
  ie lennard-jones or electrostatic or ewald etc...
*/
#include "mdmol.h++"

//#define VERBOSE
//#define NOELEC

static int mnow=0; //!< static varible which hold the static memory allocated
static int *idx=NULL; //!< pointer to interaction # 
static double *xi=NULL; //!< pointer to positions 
static double *qi=NULL; //!< pointer to charges 
static double potSRs; //!< the short range potential energy
static double potLRs; //!< the long range potential energy
static double potSRLJ; //!< the short range lj energy
static double potSREl; //!< the short range electrostatic energy.
static double potLRElk; //!< the long range electrostatic k-space energy
static double potLRElc; //!< long range electrostatic energy correction

/*! 
  routine frees the scratch memory used in energy routines
  not really needed, but good to check for memory leaks.
*/
void freeEnergyScr(void)
{ 
  if(mnow != 0){ free(xi);  free(qi); free(idx); mnow=0;}
}

/*!
  routine checks if we need more scratch
  makes sure that we have at least 
  \param mem bytes in the scratch
*/
static void checkForceScr(int mem)
{
  if(mnow < mem){
    freeEnergyScr();
    mnow = mem;
    xi = (double *)malloc(mnow*DIM*4*sizeof(double));
    qi = (double *)malloc(mnow*2*sizeof(double));
    idx = (int *)malloc(mnow*sizeof(int));
  }
}

void gatherxi(int nintr,int ip[],double qi[],double q[],int idx[],int indx[],
	      int ntypes,double xi[],double x[])
{
  int i;
  for(i=0;i<nintr;i++){
    int k;      
    int ii = ip[i*2];      
    int jj = ip[i*2+1];
    qi[i*2] = q[ii];      
    qi[i*2+1] = q[jj];
    idx[i] = indx[ii]*ntypes+indx[jj];
    ii *= DIM;
    jj *= DIM;   
    int i2 = i*2*DIM; int i2d = i2+DIM;
    for(k=0;k<DIM;k++){xi[i2+k]=x[ii+k];xi[i2d+k]=x[jj+k];}
#ifdef VERBOSE
    int m=i*DIM*2;
    printf("VERBOSE:%3d %3d %3d %4g %4g",i,ip[i*2],ip[i*2+1],qi[i*2],qi[i*2+1]);
    for(k=0;k<DIM;k++) printf("%13g %13g",xi[m+k],xi[m+DIM+k]);printf("\n");
#endif // VERBOSE
  }
}

/*! 
  calculate the short range potential energy
  \param coords is the coordinates of the systems
  \param ngbr is the neighbour lists of the systems
*/
double potSR(COORDS &coords,NEIGHBOR &ngbr)
{
  int nintr;
  int *indx = coords.itype();
  double *x = coords.x();
  double *q = coords.qch();
  int ntypes = coords.ntypes();
  int ipair = ngbr.getipair();
  int done=0;
  double poteng = 0.;
  int tinter = 0;
  ngbr.setintr();
  while(!done){
    done = ngbr.interact();
    nintr=ngbr.getintrnow();
    checkForceScr(nintr);
    int *ip = ngbr.getip();
    gatherxi(nintr,ip,qi,q,idx,indx,ntypes,xi,x);
    tinter += nintr;
    poteng += pot_lj(nintr,xi,idx);
    poteng += potElec(nintr,xi,qi);
  }
#ifdef PARA
  double rpot=poteng;
  MPI_Reduce(&rpot,&poteng,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif
#ifdef VERBOSE
  printf("Ninter = %d\n",tinter);
#endif
  
  if (ipair==1 || ipair==3 || ipair==5) poteng /= 2.;
  return poteng;
}

/*! 
  calculate the long range potential energy
  \param coords is the coordinates of the system
  \param ngbr is the neighbor call which hold the neighbor
  list information
*/
double potLR(COORDS &coords,NEIGHBOR &ngbr)
{
  double potlr = 0.;
  if(getIperd()==3){
    potlr = k_ewald(coords.natoms(),coords.x(),coords.qch());
    potlr += potEcorr(coords.x(),coords.qch());
  }
  return potlr;
}

/*! 
  calculate the short range force
  \param coords is the coordinates of the system
  \param ngbr is the neighbor cell which holds the neighbor
  list information
*/
double forceSR(COORDS &coords,NEIGHBOR &ngbr)
{
  int i,ii,nintr;
  int ntypes = coords.ntypes();
  int *indx = coords.itype();
  double *q = coords.qch();
  double *x = coords.x();
  double *a = coords.a();
  int ipair = ngbr.getipair();
  potSRLJ = potSREl = 0.;
  int done = 0;
  int tinter = 0;
  ngbr.setintr();
  while(!done){
    done = ngbr.interact();
    checkForceScr(nintr=ngbr.getintrnow());
    
    int *ip = ngbr.getip();
    gatherxi(nintr,ip,qi,q,idx,indx,ntypes,xi,x);

    tinter += nintr;
    ii = nintr*2*DIM; // copy tmp vector since force destroys it.
    for(i=0;i<ii;i++) xi[ii+i] = xi[i];
    potSRLJ += force_lj(nintr,xi,idx);
#ifndef NOELEC
    potSREl += forceElec(nintr,&(xi[ii]),qi);
    for(i=0;i<ii;i++) xi[i] += xi[ii+i]; // add both accells back to xi
#endif
    for(i=0;i<nintr;i++){
      int ii = ip[i*2]*DIM; int ix = i*DIM*2;
      for(int k=0;k<DIM;k++) a[ii+k] += xi[ix+k];
    }
    if(ipair==0 || ipair==2 || ipair==4) { // use newton's second law
      for(i=0;i<nintr;i++){
	int jj = ip[i*2+1]*DIM; int ix = i*DIM*2+DIM;
	for(int k=0;k<DIM;k++) a[jj+k] += xi[ix+k];
      }
    } // end if
  }
#ifdef PARA
  // always need to communicate if ATOMLINEAR is not defined
  if(ipair==0 || ipair==2 || ipair==4) { 
    double *xni = ngbr.getrecbuf();
    /* for(i=0;i<nmol*DIM;i++) recbuf[i] = 0.; */
    MPI_Allreduce(a,xni,DIM*coords.natoms(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<DIM*coords.natoms();i++) a[i] = xni[i];
  }
#endif

#ifdef VERBOSE
  printf("Ninter = %d\n",tinter);
#endif
  if (ipair==1 || ipair==3 || ipair==5) {potSRLJ /= 2.;potSREl /= 2.;}
  return potSRLJ+potSREl;
}

/*! calculate the long range force
  \param coords is the coordinates of the systems
  \param ngbr is the neighbour lists of the systems
  \return the total long range potential energy
*/
double forceLR(COORDS &coords,NEIGHBOR &ngbr)
{
  potLRElk=potLRElc=0;
  if(getIperd()==3){
    double *x = coords.x();
    double *a = coords.a(); 
    double *qch = coords.qch();
    potLRElk += fk_ewald(coords.natoms(),x,a,qch);
    potLRElc += ecorr(x,a,qch);
  }
  return potLRElk+potLRElc;
}

/*!
  check the forces numerically and analytically
  \param coords is the coordinates of the systems
  \param ngbr is the neighbour lists of the systems
  exits after printing out the values
*/
void force_check(COORDS &coords,NEIGHBOR &ngbr)
{
  int i,ii;
  int n = coords.natoms();
  double *x = coords.x();
  double *a = coords.a();
  
  forceSR(coords,ngbr);
  forceLR(coords,ngbr);
  printf("potSR = %g, potLR = %g\n",potSR(coords,ngbr),potLR(coords,ngbr));
  for(i=0;i<n;i++){
    int k;
    double h = 1.e-6;
    double pm,pp,xo;
    for(k=0;k<DIM;k++){
      ii = i*DIM+k;
      xo = x[ii];
      x[ii] = xo - h;  pm = potLR(coords,ngbr)+potSR(coords,ngbr);
      x[ii] = xo + h;  pp = potLR(coords,ngbr)+potSR(coords,ngbr);
      x[ii] = xo;
      printf("%4d %4d %18.10g %4d %18.10g %18.10g\n",
	     i,k,(pm-pp)/(2*h)-a[ii],ii,(pm-pp)/(2*h),a[ii]);
    }
  }
  double *qch = coords.qch();
  double d[3];
  cellDipole(n,x,qch,d);
  double cor = (4./3.)*M_PI;
  printf("pre cor = %g d=",cor);
  printf(" %g %g %g\n",d[0],d[1],d[2]);
  exit(1);
}


/*! 
  calculate the potential energy
  either by direct method or previous force calculation
  \param coords is the coordinates of the systems
  \param ngbr is the neighbour lists of the systems
  \param icalc if icalc is 1 then calculate the potentail again
  if icalc is 0 then report the previously calculated potential
*/
double pot(COORDS &coords,NEIGHBOR &ngbr,int icalc)
{
  double poteng=0;
  if(icalc==1){
    poteng += potSR(coords,ngbr); // short range non-bonded
    poteng += potLR(coords,ngbr); // long range non-bonded
  } else {
    poteng += potSRs + potLRs;
  }
  return poteng;
}

/*! 
  gets the potential energies from the static variables 
  \param sr is the short range 
  \param lr is the short range 
  \param plr is the lennard jones force
  \param pelsr is the short range electrostatic 
  \param pellr is the long range electrostatic  (ie k-space part)
  \param pecoor is the correction to the long range electrostatic 
*/
double getpot(double *sr,double *lr,double *plj,
	      double *pelsr,double *pellr,double *pecoor)
{
  *sr = potSRs;
  *lr = potLRs;
  *plj = potSRLJ;
  *pelsr = potSREl;
  *pellr = potLRElk;
  *pecoor = potLRElc;
  return potSRs + potLRs;
}

/*! 
  calculate the force on each particle
  \param coords is the coordinates of the systems
  \param ngbr is the neighbour lists of the systems
*/
void force(COORDS &coords,NEIGHBOR &ngbr)
{
  int i;
  double *a = coords.a();
  for(i=0;i<DIM*coords.natoms();i++) a[i] = 0.0;
  
  // cout<<directSum(coords.natoms(),coords.x(),coords.qch()) << endl;exit(1);
  // force_check(coords,ngbr);
  // double d[3]; d[0]=.5;d[1]= 0.;d[2]=.0; //d[1] = d[2] = .5;
  // potElecx(coords.natoms(),coords.x(),coords.qch(),d); exit(1);

  potSRs = forceSR(coords,ngbr); // short range non-bonded
  potLRs = forceLR(coords,ngbr); // long range non-bonded
  zeroFreeze(coords.a());
  zeroFreeze(coords.v());
}

