/*! 
  /file interpolation.cc
  file with routines to calculate LJ energies/forces/hessian 
  eventaully any interpolation that need to be done
*/

#include "mdmol.h++"

//#define VERBOSE

static int nxdis=0;
static int ntypes=0;
static double *xdis=NULL;
static double **finteract;

/*! 
  freeInterScr frees the scratch space used by interpolation routines
*/
void freeInterScr(void){
  if(nxdis != 0) free(xdis);
  nxdis = 0;
  // if(xdis != NULL) delete[] xdis;
}

/*! 
  checkInterScr check that there is enough scratch space and if
  not then adds some to the scratch space
*/
static void checkInterSrc(int nintr){
  if(nintr>nxdis){
    freeInterScr();
    xdis = (double *)malloc(nintr*DIM*sizeof(double));
    // xdis = new double[nintr*DIM];
    nxdis = nintr;
  }  
}

static void LJ_pot(double r,double eps,double sig, double func[3])
{
  double sor6 = r/sig;
  sor6 *= sor6;
  sor6  = sor6*sor6*sor6;
  func[0] = 4.0*eps*sor6*(sor6-1.);
  func[1] = 24.0*eps*sor6*(1.-2.*sor6); // virial
  func[2] = 24.0*eps*sor6*(26.*sor6 - 7.)/(r*r); // second derivative
  
  return pot;
}

/*! 
   InterSetup sets up Interpolation interactions
   \param ntypes is the number of different typs
   \param eps is a double array [ntypes] which are the different LJ well depths
   \param sig is a double array [ntypes] which are the different LJ Lenths
   \param rcut is a double with cut for the LJ interactions
   \param xbox is the box lenth
   cross terms are calculated with the combining rules, <br>
   eps[i][j] = sqrt(eps[i]*eps[j]), sig[i][j] = .5*(sig[i]+sig[j])
*/
void Intersetup(int ntyps,double eps[],double sig[],double rcut,double xbox)
{
  // eps = 0; no LJ
  xbox=1.; // don't scale LJ sigma's interactions
  ntypes = ntyps;
  int i,j;
  for(i=0;i<ntypes;i++) sig[i] /= xbox; // scale distance to xbox
  // allocate memory
  int ntypes2 = ntypes*ntypes;
  eps4 = (double *)malloc(2*ntypes2*sizeof(double));
  sig6 = eps4+ntypes2;
  // set up interactions
  for(i=0;i<ntypes;i++){
    for(j=0;j<ntypes;j++){
      eps4[i*ntypes+j] = 4.*sqrt(eps[i]*eps[j]);
      double sigij = .5*(sig[i]+sig[j]);
      sig6[i*ntypes+j] = sigij*sigij*sigij*sigij*sigij*sigij;
    }
  }
  // sig6[0] = 1.5625e-8;
  rcut2 = rcut*rcut/(xbox*xbox);
}

/*! 
  pot_lj calculates the LJ potential energy between particles
  \param nintr is of interactions
  \param xi is a double array [nintr*2*DIM] which has the particle positions
  that are interacting the format is x0 x1 x2 x3 x4 x5 .... x[nintr*2*DIM] 
  where x0 x1 x2 is the x, y, z position of atom 1 which is interaction with
  x3, x4, x5 (x,y,z) positions of atom 2 etc...
  \param idx is an integer array [nintr] which hold the interaction type
  for the ith interaction
  the returns the total potential energy
*/
double pot_inter(int nintr,double xi[],int idx[])
{
  int i,ii,k;
  double pot=0.;

  checkLJSrc(nintr);

  for(i=0;i<nintr;i++){
    ii = i*DIM*2;
    for(k=0;k<DIM;k++) xdis[k+i*DIM] = xi[ii+k] - xi[ii+DIM+k];
  }
  periodic(nintr,xdis);
  for(i=0;i<nintr;i++){
    double sor6,dis2=0.;
    for(k=0;k<DIM;k++) {sor6 = xdis[k+i*DIM];dis2 += sor6*sor6;}
    if(dis2<rcut2){
      switch(type[idx[i]]){
      case 0:  // null interaction 
	break;
      case 1: // LJ
	k = idx[i];
	sor6  = sig6[k]/(dis2*dis2*dis2);
	pot  += eps4[k]*sor6*(sor6-1.);
	break;
	// virial += eps24*sor6*(1.-2.*sor6);
      case 2: // Wilson
	k = idx[i]; 
	// A*exp(-a (x-x0)) - sig/r
	pot += eps4[k]*sor6*sor6;
      case 3: // spline
	break;
      case 4: // 1/r
	break;
      default:
	error();
	break;
      }
    }
  }
  return pot;
}

/*--------------------------------------------------------*/
/*! force_lj
  calculates the force and potential energy between particles 
  \param nintr is an int for the number of interactions
  \param xi is a double array [nintr*2*DIM] which has the particle positions
  that are interacting the format is x0 x1 x2 x3 x4 x5 .... x[nintr*2*DIM] 
  where x0 x1 x2 is the x, y, z position of atom 1 which is interaction with
  x3, x4, x5 (x,y,z) positions of atom 2 etc...
  \param idx is an integer array [nintr] which holds the interaction type
  for the ith interaction.  
  on exit force_lj returns the total potential energy, and the forces
  of the interactions in the xi array. (ie the xi array is destroyed and
  replaced with the forces)
*/
double force_inter(int nintr,double xi[],int idx[])
{
  int i,ii,id,k;

  checkLJSrc(nintr);
  double pot = 0;

  for(i=0;i<nintr;i++){
    ii = i*DIM*2;
    for(k=0;k<DIM;k++) xdis[k+i*DIM] = xi[ii+k] - xi[ii+DIM+k];
  }
  periodic(nintr,xdis);
  for(i=0;i<nintr;i++){
    id = i*DIM;
    ii = id*2;
    double dis2=0.0,sor6,virial,fcc;
    for(k=0;k<DIM;k++) {double dx = xdis[k+id];dis2 += dx*dx;}
    if(dis2<rcut2){
      k = idx[i];
      sor6  = sig6[k]/(dis2*dis2*dis2);
#ifdef REPULSIVE
      virial = -12.*eps4[k]*sor6*sor6; 
      pot += eps4[k]*sor6*sor6;
#else // LJ
      virial= 6.*eps4[k]*sor6*(1.-2.*sor6);
      pot += eps4[k]*sor6*(sor6-1.);
#endif
      // sprs += virial;
      fcc  = virial/dis2;      
      for(k=0;k<DIM;k++) {
	double frc = xdis[k+id]*fcc;
	xi[ii+k] = -frc; 
	xi[ii+DIM+k] = frc;
      }
    } else for(k=0;k<DIM*2;k++) xi[ii+k] = 0.0; // zero acceleration
  }
  // double rho=nmol/xbox;
  // *tpot = *tpot+vlrc0*rho;
  // *sprs = *sprs/3.+wlrc0*rho;
  return pot;
}
/*--------------------------------------------------------*/
/*! d2_lj
  calculates the hessian (force constant matrix) from the LJ interaction
  \param natoms an int of number of particles
  \param x is a double array [natoms*DIM] which has the particle positions
  \param type is int array [natoms] which has the atom types
  \param fmat is a double array [natoms][natoms]
  d2_lj fills the force matrix fmat
*/
void d2_lj(int nintr,double x[],int type[],NEIGHBOR &ngbr,double **fmat)
{
  static int istart=1;
  static double **fbl;
  int i,i2,ii,jj,k,l;

  // rcut2 = rcut*rcut;
  double *xi = ngbr.getxi();
  int *ip = ngbr.getip();
  if(istart){
    istart=0;
    fbl = dmatrix(0,DIM-1,0,DIM-1);
  }
  double *xdis = new double[nintr*DIM];
  for(i=0;i<nintr;i++){
    ii = i*DIM*2;
    for(k=0;k<DIM;k++) xdis[k+i*DIM] = xi[ii+k] - xi[ii+DIM+k];
  }
  periodic(nintr,xdis);
  for(i=0;i<nintr;i++){
    double du,d2u,duor,dfac;
    double dis2 = 0.;
    for(k=0;k<DIM;k++) {double dx = xdis[k+i*DIM];dis2 += dx*dx;}
    double r = sqrt(dis2);
    if(dis2>rcut2){
      du = d2u = 0.;
    } else {
      k = type[ip[i*2]]*ntypes+type[ip[i*2+1]];
      double sor6 = sig6[k]/(dis2*dis2*dis2);
      du   = 6.*eps4[k]*sor6*(1.-2.*sor6)/r;
      d2u  = 6.*eps4[k]*sor6*(26.*sor6 - 7.)/dis2;
    }
    duor = du/r;
    dfac = duor - d2u;
    for(k=0;k<DIM;k++) xdis[k+i*DIM] /= r;
    
    for(k=0;k<DIM;k++) {
      for(l=0;l<DIM;l++) fbl[k][l] = dfac*(xdis[k+i*DIM]*xdis[l+i*DIM]);
      fbl[k][k] -= duor;
    }
    i2 = i*2;
    ii = ip[i2];
    jj = ip[i2+1];
    for(k=0;k<DIM;k++) {
      for(l=0;l<DIM;l++) {
	fmat[jj+l][ii+k] = fmat[ii+k][jj+l] = fbl[k][l];
	/* add contributions to diagonal block */
	fmat[ii+k][ii+l] -= fbl[k][l];
	fmat[jj+k][jj+l] -= fbl[k][l];
      }  
    }
  }
  delete[] xdis;
}
