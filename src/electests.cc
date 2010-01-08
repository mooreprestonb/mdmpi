/*! \file electests.cc
   subroutines included in the charge-charge interactions 
   test for real space and ewald sums
*/

/*!
  old test routine to calculate the force and virial between to 
  all the pairs of atoms real space terms using the ewald summation.
*/
void fr_ewald(int nch,double *px,double *fx,double *qch)
{
  int i,j,k,ii,jj;
  double qi,dx[3],fxi[3],r,r2,sqrt_pi;
  double expkr2,erfkr,vir;

  sqrt_pi = sqrt(M_PI);
  
  for(i=0;i<nch-1;i++){
    ii = i*DIM;
    qi = qch[i];
    for(k=0;k<DIM;k++) fxi[k] = 0.;
    for(j=i+1;j<nch;j++){
      jj = j*DIM;
      for(k=0;k<DIM;k++) dx[k] = px[ii+k]-px[jj+k];
      periodic(1,dx);

      r2 = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
      r = sqrt(r2);
      erfkr = erfc(kappa*r);
      expkr2 = exp(kappa2*r2);
      vir = qi*qch[j]*(erfkr/r+2.*kappa/(sqrt_pi*expkr2));
      for(k=0;k<DIM;k++) dx[k] *= vir/r2;
      for(k=0;k<DIM;k++) fxi[k] -= dx[k];
      for(k=0;k<DIM;k++) fx[jj+k] += dx[k];
    }
    for(k=0;k<DIM;k++) fx[ii+k] += fxi[k];
  }
}
/*-------------------------------------------------------------------*/
/*! calculate the 3d periodic potential using direct summation */
double directSum(const int nch,const double px[],const double qch[])
{
  int i,j,k,l,m,n;
  double spot;
  double R[3];
  double box = getBox();
  int iperd = getIperd();
  setIbox(0,box);
  double *p1x = new double[nch*DIM];
  double *qi = new double[2*nch];
  double *p3x = new double[2*nch*DIM];
  
  for(i=0;i<nch*DIM;i++) p1x[i] = px[i];
  periodic(nch,p1x);
  double rcut2s = rcut2;
  int kmax = 40;
  rcut2 = 3.*(kmax+1)*(kmax+1)*box*box;
  
  for(j=0;j<nch;j++) qi[2*j+1] = qch[j];
  spot = 0.;  
  for(l=-kmax;l<=kmax;l++){
    R[0] = l*box;
    for(m=-kmax;m<=kmax;m++){
      R[1] = m*box;
      for(n=-kmax;n<=kmax;n++){
	R[2] = n*box;
	if(n==0 && m==0 and l==0) continue;
	for(j=0;j<nch;j++){	  
	  int jj = j*DIM;
	  for(k=0;k<DIM;k++) p3x[jj*2+k+DIM] = p1x[jj+k]+R[k]; 
	}
	for(i=0;i<nch;i++){
	  for(j=0;j<nch;j++) qi[2*j] = qch[i];
	  int jj = i*DIM;
	  for(j=0;j<nch;j++) for(k=0;k<DIM;k++) p3x[j*2*DIM+k] = p1x[jj+k];
	  if(n==0 && m==0 and l==0) {
	    qi[i*2] = 0;  // don't compute self
	    p3x[i*DIM*2] += 1; // don't divide by zero
	  }
	  spot += potElec(nch,p3x,qi);
	}
      }
    }
  }
  rcut2 = rcut2s;
  setIbox(iperd,box);
  delete[] p1x;
  delete[] p3x;
  delete[] qi;
  return .5*spot;
}

/*! 
  calculate the potential at a particular point 
*/
void potElecx(int nch,double *xi,double *qi, double x[])
{
  int i,ii,k;
  double vr=0.;
 
  checkEwaldScr(nch);
  
  for(i=0;i<nch;i++){
    ii = i*DIM;
    for(k=0;k<DIM;k++) xdis[k+i*DIM] = xi[ii+k] - x[k];
  }
  // periodic(nch,xdis);
  for(i=0;i<nch;i++){
    double dx,dis2=0.;
    for(k=0;k<DIM;k++) {dx = xdis[k+i*DIM];dis2 += dx*dx;}
    // if(dis2<rcut2){
      double r=sqrt(dis2);
      printf("%d %g %g\n",i,qi[i],r);
      if(r==0) fprintf(stderr,"ERROR zero distance\n");
      else vr += qi[i]/r;
      // }
  }
  printf("pot at z = %g\n",vr);
  
  double vol=1;
  double d[3];
  cellDipole(nch,xi,qi,d);
  printf("dipole correction = %g\n",
	 (4.*M_PI/(3.*vol))*(x[0]*d[0]+x[1]*d[1]+x[2]*d[2]));
  
  double q = cellQuadrupole(nch,xi,qi);
  printf("quad correction = %g\n",2*M_PI*q/3.);
  exit(1);
}
