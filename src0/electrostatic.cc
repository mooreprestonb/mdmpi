/*! \file electrostatic.cc
   subroutines included in the charge-charge interactions 
   (ie 1/r and ewald sums) 
   see for example TM Nymand and P Linse JCP V112 p6152 (2000)
*/

/* 
   Copyright (C) 1997-2003 Dr. Preston B. Moore

Dr. Preston B. Moore
Associate Director, Center for Molecular Modeling (CMM)
University of Pennsylvania, Department of Chemistry, Box 188 
231 S. 34th St. Philadelphia, PA 19104-6323 USA
EMAIL: moore@cmm.chem.upenn.edu  
WWW: http://www.cmm.upenn.edu/~moore

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include "mdmol.h++"

/* #define CUBIC */
/* #define ANALYTIC */
/* #define CHECK_BKG */
/* #define PCORR */
/* #define DEBUG */

/* #define DEBUG */
// #define RCUTOFF
#define TPI  (2.*M_PI)
#define HELPFULL

#ifdef ANALYTIC
extern double alpha_ewald;
#endif

static int icutoff;
static int ktot,**kv=NULL;
static double kappa,kappa2,falp2,self,back;
static double *cosscr,*sinscr,*helrscr,*heliscr;
static double vk_ewald=0,vp_ewald=0,vp_ewald_tensor[9];
static double pot_rs=0;
static int ncorr=0,*icorr=NULL;
static double pot_ecorr=0,prs_ecorr=0,prs_ecorr_tensor[9];
static int nxdis=0;
static double *xdis=NULL;
static double rcut2;
static double esuf=-1;

/*! 
  return ewald electrostatic potential 
  which has been stored while calculation the force
*/
double eleclrpot(void){return vk_ewald;} 
/*! 
  return real space electrostatic potential 
  which has been stored while calculating the force
*/
double elecsrpot(void){return pot_rs;} 
/*! 
  return ewald self correction electrostatic potential 
  which has been stored while calculating the force
*/
double eleccorrpot(void){return pot_ecorr;}

/*! 
  free the ewald scratch space 
*/
void freeEwaldScr(void) { if(xdis!=NULL) {nxdis=0;free(xdis);}}

/*! 
  check to see if we have enough scratch space and 
  if not allocate some more
*/
static void checkEwaldScr(int nintr){
  if(nintr>nxdis){
    freeEwaldScr();
    xdis = (double *)cmalloc(nintr*DIM*sizeof(double));
    // xdis = new double[nintr*DIM];
    nxdis = nintr;
  } 
}

/*! 
   free the static variable and scratch space used in 
   the electrostatic calculations
*/
void freeEwaldStatic(void)
{
  if(kv != NULL){
    int i;
    for(i=0;i<DIM+1;i++) free(kv[i]);
    free(kv); free(cosscr);
  }
  freeEwaldScr();
  if(icorr!=NULL){ ncorr=0;free(icorr);}
}

/*!
  ewald setup, takes the charges some parameters and sets up the 
  k vectors etc...
  must be called before ewald is called
*/
void ewaldSetup(int nch,double *qch,int kmax,double alpha,double rcut,
		double esufs,int perd,int cutoff)
{
  int i,j,k;
  int ksqmax,memk,ksq,membytes=0;
  double sqrt_pi;
  char line[MAXLINE];

  // static members

  rcut2 = rcut*rcut;
  sqrt_pi = sqrt(M_PI);
  kappa = alpha;
  kappa2 = kappa*kappa;
  falp2 = 1./(4.*kappa2);
  int iperd = perd;
  icutoff = cutoff;
  esuf = esufs;

  self = 0.;
  back = 0.;

  for(j=0;j<nch;j++){
    self += qch[j]*qch[j];
    back += qch[j];
  }
  if(iperd != 3 || icutoff==1 || self == 0){
    ktot = 0;
  } else {
    self *= -kappa/sqrt_pi;
    back = -.5*back*back*M_PI/kappa2;
#ifdef CHECK_BKG
    if(back < ERRMAX) back = 0.;
#endif
    memk = kmax*kmax*kmax*kmax;
    kv = (int **)cmalloc((DIM+1)*sizeof(int *));
    for(k=0;k<DIM;k++) kv[k] = (int *)cmalloc(memk*sizeof(int));
    kv[DIM] = (int *)cmalloc((memk+1)*sizeof(int));
    
    ktot = 0;
    ksqmax = kmax*kmax;
    
    // sprintf(line,"Setting up k-vectors for ewald sum (kmax = %d)",kmax);
    // writeOut(line);
    for(i=0;i<=kmax;i++){
      for(j = (i==0)?0:-kmax;j<=kmax;j++){
        kv[DIM][ktot] = 1;
        for(k= (i==0 && j==0)?1:-kmax;k<=kmax;k++){
          ksq = i*i+j*j+k*k;
          if(ksq<=ksqmax) {
            kv[0][ktot] = i; kv[1][ktot] = j; kv[2][ktot] = k;
            kv[DIM][ktot+1] = 0;
            ktot++;
            if(memk==ktot){
              memk += kmax*kmax*kmax;
	      for(k=0;k<DIM;k++)
		kv[k] = (int *)realloc(kv[k],memk*sizeof(int));
              kv[DIM] = (int *)realloc(kv[DIM],(memk+1)*sizeof(int));
            }
          }
        }
      }
    }
    
    for(k=0;k<DIM;k++) kv[k] = (int *)realloc(kv[k],ktot*sizeof(int));
    kv[DIM] = (int *)realloc(kv[DIM],(ktot+1)*sizeof(int));
    membytes += (4*ktot+1)*sizeof(int);
    cosscr = (double *)cmalloc(4*nch*sizeof(double));
    sinscr = cosscr +   nch;
    helrscr  = cosscr + 2*nch;
    heliscr  = cosscr + 3*nch;
    membytes += 4*nch*sizeof(double);
  }

  sprintf(line,"There were %d k-vectors set up.",ktot);writeOut(line);
  sprintf(line,"Self correction term = %g K",self); writeOut(line);
  sprintf(line,"Back ground term = %g K",back); writeOut(line);
  sprintf(line,"Bytes of memory allocated for electrostatic = %d",membytes);
  writeOut(line);
}

/*!
  calculate the cell dipole
*/
void cellDipole(int nch,double *px,double *qch,double d[DIM])
{
  int j,k;
  double dx[DIM];
  for(k=0;k<DIM;k++) d[k] = 0;
  for(j=0;j<nch;j++){
    double qi = qch[j]; int jj = j*DIM;
    for(k=0;k<DIM;k++) dx[k] = px[jj+k]; periodic(1,dx);
    for(k=0;k<DIM;k++) d[k] += qi*dx[k];
  }
}

/*!
  calculate the cell quadrupole
*/
double cellQuadrupole(int nch,double *px,double *qch)
{
  int i,k;
  double Q=0;
  double dx[DIM];
  for(i=0;i<nch;i++){
    double tr=0; int l = i*DIM;
    for(k=0;k<DIM;k++) dx[k] = px[l+k]; periodic(1,dx);
    for(k=0;k<DIM;k++) tr += dx[k]*dx[k];
    Q += qch[i]*tr;
  }
  return Q;
}

/*! 
  calculate the surface term for the ewald sums
  use esuf==-1 as esuf==infinity
*/
double surf_ewald(int nch,double *xi,double *qch)
{
  if(esuf==-1) return 0;
  int k;
  double pref = 2.*M_PI/((2.*esuf+1.)*getVol());
  double suf=0;
  double du[3];
  cellDipole(nch,xi,qch,du);
  for(k=0;k<DIM;k++) suf += du[k]*du[k];
  return pref*suf;
}

/*! calculate the surface force for the ewald sums
  use esuf==-1 as esuf==infinity
*/
void surf_fewald(int nch,double *xi,double *fx,double *qch)
{
  if(esuf==-1) return;
  int j,k;
  double pref = -4.*M_PI/((2.*esuf+1.)*getVol());
  double du[3];
  cellDipole(nch,xi,qch,du);
  for(k=0;k<DIM;k++) du[k] *= pref;
  for(j=0;j<nch;j++){
    double qi = qch[j]; int jj = j*DIM;
    for(k=0;k<DIM;k++) fx[jj+k] += qi*du[k];
  }
}
/*!
  calculate the potential for the k-space potential 
*/
double k_ewald(int nch,double *px,double *qch)
{
  int i,j,k,jj;
  double ak[3],vk,qi;
  double sumr,sumi,g2,preg,smag,rk;
  double hmat[9],hmati[9];

  double xbox = getBox();
  for(i=0;i<9;i++) hmati[i] = hmat[i]=0;
  hmat[0] = hmat[4] = hmat[8] = xbox;
  hmati[0] = hmati[4] = hmati[8] = 1./xbox;

  double pref = 2.*2.*M_PI/getVol(); // mult by 2 since we have only positive 
  vk = 0.;                           // half of k vect.

  for(i=0;i<ktot;i++){
    for(k=0;k<DIM;k++) 
      ak[k]=(kv[0][i]*hmati[0+k]+kv[1][i]*hmati[3+k]+kv[2][i]*hmati[6+k])*TPI; 

    g2 = ak[0]*ak[0]+ak[1]*ak[1]+ak[2]*ak[2]; // m^2*4*pi^2
    preg = pref*exp(-g2*falp2)/(g2);
    sumi = sumr = 0.;
    for(j=0;j<nch;j++){
      qi = qch[j];
      jj = j*DIM;
      for(rk=0.,k=0;k<DIM;k++) rk += px[jj+k]*ak[k];
      sumr += qi*cos(rk);
      sumi += qi*sin(rk);
    }
    smag = sumr*sumr+sumi*sumi;
    vk += preg*smag;
  }
  vk_ewald = vk+self+back/getVol() + surf_ewald(nch,px,qch);
  return vk_ewald;
}

/*!
  calculate the force and potential for the k space term
  in the ewald sum, returns the potential as well as storing
  the potential for other calls.
*/
double fk_ewald(int nch,double *px,double *fx,double *qch)
{
  int i,j,k;
  double rvol,vk,vp,q;
  double sumr,sumi,g2,preg,prep,smag,rk;
  double sr[3],si[3],ak[3],hmati[9],hmat[9];
  double vk_local,vp_local,vp_local_tensor[9],vp_tensor[9];
  double rhr,rhi;
#ifdef PARA
    // xdis = new double[nintr*DIM];
  double *scr_buf = new double[2*nch*DIM];
  double *scr_rec = scr_buf + nch*DIM;    
#endif

  double xbox = getBox();
  rvol = 1./getVol();
  vp = vk = vp_local = vk_local=0.;
  for(i=0;i<9;i++) vp_tensor[i] = vp_local_tensor[i]=0.;

#ifdef PARA
  for(i=0;i<nch*DIM;i++) scr_buf[i]=scr_rec[i]=0.0;
#endif

  for(i=0;i<9;i++) hmati[i] = hmat[i]=0;
  hmat[0] = hmat[4] = hmat[8] = xbox;
  hmati[0] = hmati[4] = hmati[8] = 1./xbox;

  /* get sc compenents for recursion */

  if(ktot != 0){
    for(j=0;j<nch;j++){
      int jj = DIM*j;
      rk = (px[jj]*hmati[6] + px[jj+1]*hmati[7] + px[jj+2]*hmati[8])*TPI;
      cosscr[j] = cos(rk); sinscr[j] = sin(rk);
    }
  }

  for(i=0;i<ktot;i++){
    si[0] = TPI*kv[0][i];    si[1] = TPI*kv[1][i];    si[2] = TPI*kv[2][i];
    for(k=0;k<DIM;k++)ak[k]=(si[0]*hmati[k]+si[1]*hmati[k+3]+si[2]*hmati[k+6]);

    g2 = ak[0]*ak[0]+ak[1]*ak[1]+ak[2]*ak[2];
    preg = 2.*exp(-g2*falp2)*rvol/g2;
    /* preg = exp(-g2*falp2)/g2; */
    prep = -2.*preg*(g2*falp2+1.)/g2;

#ifndef HELPFULL
    kv[3][i] = 1; test to see that recursion works
#endif
    if(kv[3][i] || i==0){
      for(j=0;j<nch;j++){
	int jj = j*DIM;
	rk = (px[jj]*ak[0] + px[jj+1]*ak[1] + px[jj+2]*ak[2]);
	helrscr[j] = cos(rk);  heliscr[j] = sin(rk);
      }
    } else {
      /* update helpfull vectors with recursion 
         ie cos(a+b) = cos(a)cos(b)-sin(a)sin(b) etc..*/
      for(j=0;j<nch;j++){
        rhr = helrscr[j];   rhi = heliscr[j];
        sumr = cosscr[j]; sumi = sinscr[j];
        helrscr[j] = rhr*sumr-rhi*sumi;
        heliscr[j] = rhi*sumr+rhr*sumi;
      }
    }
    sumr = sumi = 0.;
    for(j=0;j<nch;j++){
      q = qch[j]; 
      sumr += q*helrscr[j]; 
      sumi += q*heliscr[j];
    }
    smag = sumr*sumr+sumi*sumi;

    vk_local += 2.*M_PI*preg*smag;
    vp_local += prep*smag*g2;
    
    for(k=0;k<DIM;k++){
      for(j=0;j<DIM;j++){
	vp_local_tensor[j+k*DIM] += prep*smag*ak[k]*ak[j];
      }
    }
    sumr *= 4.*M_PI*preg; 
    sumi *= 4.*M_PI*preg;

    for(j=0;j<nch;j++){
      q = qch[j]; 
      double qr = q*sumr;
      double qi = q*sumi;
      rhr = helrscr[j]; rhi = heliscr[j];
      for(k=0;k<DIM;k++) {sr[k] = ak[k]*qr; si[k] = ak[k]*qi;}
#ifdef PARA
      for(k=0;k<DIM;k++) scr_buf[DIM*j+k] +=  sr[k]*rhi - si[k]*rhr;
#else
      for(k=0;k<DIM;k++) fx[DIM*j+k] +=  sr[k]*rhi - si[k]*rhr;
#endif
    }
  }

  surf_fewald(nch,px,fx,qch);

#ifdef PARA
  MPI_Allreduce(&vk_local,&vk,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&vp_local,&vp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(vp_local_tensor,vp_tensor,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else 
  vk = vk_local;
  vp = vp_local;
  for (i=0;i<9;i++) vp_tensor[i]=vp_local_tensor[i];
#endif

  rk = back*rvol;
  vk_ewald = vk+rk+self+surf_ewald(nch,px,qch);
  vp_ewald = vp + 3.*(vk+rk);
  
  for (i=0;i<9;i++) vp_ewald_tensor[i]=vp_tensor[i];
  vp_ewald_tensor[0] += vk+rk; // add diagional contribution
  vp_ewald_tensor[4] += vk+rk;
  vp_ewald_tensor[8] += vk+rk;
 
#ifdef PARA
  MPI_Allreduce(scr_buf,scr_rec,nch*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(j=0;j<nch*3;j++) fx[j] += scr_rec[j];

  delete[] scr_buf;
#endif
  return vk_ewald;
}

/*!
  gets the diffent contributions to the potential and pressure of
  the ewald summation
*/
void getvewald(double *pot_ewald,double *prs_ewald,double prs_ewald_tensor[9])
{
  int i;
  *pot_ewald += vk_ewald+pot_ecorr; 
  *prs_ewald += vp_ewald+prs_ecorr;
  for(i=0;i<9;i++) {
    prs_ewald_tensor[i] += (vp_ewald_tensor[i]+prs_ecorr_tensor[i]); 
  }
#ifdef PCORR
  printf("ewald vk = %lg, ecorr = %lg\n",vk_ewald,pot_ecorr);
  printf("ewald vp = %lg, ecorr = %lg\n",vp_ewald,prs_ecorr);
#endif
}

/*!
  add to particles to the ecorr part.  ie this particles do NOT
  have and interaction in within the same box
*/
void addEcorr(int i,int j)
{
  int k;
  if(j<i){k=i;i=j;j=k;}
  if(icorr==NULL){ // first case
    ncorr = 1;
    icorr = (int *)cmalloc(2*sizeof(int));
  } else {
    for(k=0;k<ncorr;k++) if(icorr[k*2]==i && icorr[k*2+1]==j) return;
    ncorr++;
    icorr = (int *)realloc(icorr,ncorr*2*sizeof(int));
  }
  k = ncorr*2-2;
  icorr[k] = i;icorr[k+1] = j;
}

/*! 
   print the atom that are in the ewald correction term
*/
void printEcorr(void)
{
  printf("# Ecorr = %d\n",ncorr);
  for(int i=0;i<ncorr;i++) 
    printf("Ecorr[%d] %d-%d\n",i,icorr[i*2],icorr[i*2+1]);
}

/*! 
  routine to calculate ewald corrections to pairs of charged
  particles that don't interact via a coulomb term in the first image 
*/
double ecorr(double *px,double *fx,double *qch)
{
  int i,k,ii,jj;
  double dvecorr,sprs,spot;
  double xdis[3],dis2,dis;
  double erfr,qijr;
  double sprs_tensor[9];

  spot = sprs = 0.;
  for(i=0;i<9;i++) sprs_tensor[i]=0.0;

  double sqrt_pi = 1./sqrt(M_PI);
  
  for(i=0;i<ncorr;i++){
    ii = icorr[i*2]; jj = icorr[i*2+1];
    qijr = qch[ii]*qch[jj];
    ii *= DIM;     jj *= DIM;
    for(k=0;k<DIM;k++) xdis[k] = px[ii+k] - px[jj+k]; 
    dis2 = xdis[0]*xdis[0]+xdis[1]*xdis[1]+xdis[2]*xdis[2];
    dis = sqrt(dis2);
    qijr /= dis;

    erfr = erf(dis*kappa);
    spot += -qijr*erfr;
    dvecorr = -qijr*(-2.*kappa*exp(-kappa2*dis2)*sqrt_pi+erfr/dis);
    sprs += dis*dvecorr;

    dvecorr /= dis;

    for(k=0;k<DIM;k++){
      for(int l=0;l<DIM;l++){
	sprs_tensor[k+l*DIM]+= dvecorr*xdis[l]*xdis[k];
      }
    }
    for(k=0;k<DIM;k++) xdis[k] *= dvecorr;
    for(k=0;k<DIM;k++) fx[ii+k] += xdis[k];
    for(k=0;k<DIM;k++) fx[jj+k] -= xdis[k];
  }
  pot_ecorr = spot;
  prs_ecorr = sprs;
  for(i=0;i<9;i++) prs_ecorr_tensor[i] = sprs_tensor[i];
  return pot_ecorr;
}

/*! 
  routine to calculate ewald potential corrections to pairs of charged
  particles that don't interact via a coulomb term in the first image 
*/
double potEcorr(double *px,double *qch)
{
  int i,k,ii,jj;
  double spot,xdis[3],dis2,dis,erfr,qijr;

  spot = 0.;  
  for(i=0;i<ncorr;i++){
    ii = icorr[i*2]; jj = icorr[i*2+1];
    qijr = qch[ii]*qch[jj];
    ii *= DIM;     jj *= DIM;
    for(k=0;k<DIM;k++) xdis[k] = px[ii+k] - px[jj+k]; 
    dis2 = xdis[0]*xdis[0]+xdis[1]*xdis[1]+xdis[2]*xdis[2];
    dis = sqrt(dis2);
    qijr /= dis;

    erfr = erf(dis*kappa);
    spot += -qijr*erfr;
  }
  pot_ecorr = spot;
  return pot_ecorr;
}

/*!
  calculate the potential between two particles
  either qi qj / r or ewald screen type of qi qj erfc(k r) / r
*/
double potElec(int nch,double *xi,double *qi)
{
  int i,ii,k;
  double vr=0.;
  int iperd = getIperd();
 
  checkEwaldScr(nch);
  
  for(i=0;i<nch;i++){
    ii = i*DIM*2;
    for(k=0;k<DIM;k++) xdis[k+i*DIM] = xi[ii+k] - xi[ii+DIM+k];
  }
  periodic(nch,xdis);
  for(i=0;i<nch;i++){
    double dx,dis2=0.;
    for(k=0;k<DIM;k++) {dx = xdis[k+i*DIM];dis2 += dx*dx;}
    if(dis2<rcut2){
      double r=sqrt(dis2);
      if(iperd==3 && icutoff==0) vr+=qi[i*2]*qi[i*2+1]*erfc(kappa*r)/r;
      else vr += qi[i*2]*qi[i*2+1]/r;
    }
  }
  pot_rs = vr;
  return vr;
}

/*!
  calculate the potential and force between two particles
  either qi qj / r or ewald screen type of qi qj erfc(k r) / r
*/
double forceElec(int nch,double *xi,double *qi)
{
  int i,ii,k;
  double vr=0.;
#ifndef RCUTOFF
  double sqrt_pi = sqrt(M_PI);
#endif
  int iperd = getIperd();
  checkEwaldScr(nch);
  
  for(i=0;i<nch;i++){
    ii = i*DIM*2;
    for(k=0;k<DIM;k++) xdis[k+i*DIM] = xi[ii+k] - xi[ii+DIM+k];
  }
  periodic(nch,xdis);
  for(i=0;i<nch;i++){
    double fcc,dx,dis2=0.;
    for(k=0;k<DIM;k++) {dx = xdis[k+i*DIM];dis2 += dx*dx;}
    if(dis2<rcut2){
      double qij = qi[i*2]*qi[i*2+1];
      double r = sqrt(dis2);
      if(iperd == 3 && icutoff==0) {
	double erfkr = erfc(kappa*r);
	double expkr2 = exp(kappa2*dis2);
	fcc = qij*(erfkr/r+2.*kappa/(sqrt_pi*expkr2))/dis2;
	vr += qij*erfkr/r;
      } else {
	fcc = qi[i*2]*qi[i*2+1]/(dis2*r);
	vr += qij/r;
      }
      ii = i*DIM*2; for(k=0;k<DIM;k++) xi[ii+k] = fcc*xdis[k+i*DIM];
      ii += DIM;    for(k=0;k<DIM;k++) xi[ii+k] = -fcc*xdis[k+i*DIM];
    } else { // zero non interactions 
      ii=i*DIM*2; for(k=0;k<2*DIM;k++) xi[ii+k] =0.; 
    }
  }
  pot_rs = vr;
  return vr;
}

/*-------------------------------------------------------------------*/

