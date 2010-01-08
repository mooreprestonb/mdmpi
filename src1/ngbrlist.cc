/*! \file ngbrlist.cc
file contains neighbor list members
*/

#include "mdmol.h++"

//! use linear instead of upper triangular method for getting neighbors
#define ATOMLINEAR
// #define skin  .01
// #define rcut2 25 5^2 */

/*! get the next neighbors in the list using the verlet list O(N) */
int NEIGHBOR :: vl(void)
{
  int i,i2,ii,jj,ni;
  
  if(nintr==-2) {fprintf(stderr,"why are we here vList\n");return 1;}
  if(nintr == -1) nintr=0;
  int nih = MIN(nintr+nix,nvl); 
  for(i=nintr;i<nih;i++){
    i2 = i*2;        ii = ipvl[i2];    jj = ipvl[i2+1];
    ni = (i-nintr)*2;  ip[ni] = ii;    ip[ni+1] = jj;
  }
  nnow = nih-nintr;
  nintr = nih;
  if(nintr==nvl) {nintr = -2;return 1;}
  return 0;
}

/*! do a N^2 search for neighbors 
 \warning we should do a link-list O(N) search for neighbors
*/
void NEIGHBOR :: getVlist(double x[])
{
  int i,ii,j,nj,k,kk,ibegin,iend;
  double skin2;

  skin2 = (rcut+skin)*(rcut+skin);
  nvl = 0; // offset of interacting particles in list
#ifdef ATOMLINEAR
  decomp1d(nmol,size,rank,&ibegin,&iend);
#else
  decomp_trig(nmol,size,rank,&ibegin,&iend);
#endif

  double *xdis = new double[nmol*DIM];

  for(i=ibegin;i<iend;i++){
    ii = i*DIM;
#ifdef ATOMLINEAR
    nj = nmol;
    if(ipair==2){
      nj=nmol/2+1;
      if((nmol%2)==0 && (i>=nmol/2)) nj = nmol/2;
    }
    j = 1;
#else  //  good old i<j
    j = ((ipair==2)?i+1:0); 
    nj = nmol;
#endif
    for(kk=0;j<nj;j++){
#ifdef ATOMLINEAR
      int jj = ((j+i)%nmol);
#else 
      int jj = j;
#endif
      if(!exclude(i,jj)) {
	jj *= DIM;
	for(k=0;k<DIM;k++) xdis[k+kk*DIM] = x[ii+k] - x[jj+k];
	kk++;
      }
    }
    periodic(kk,xdis);
    
    if(nvl > memvl-nmol){
      memvl += nmol; // allocate more mem and store how much we have
      ipvl = (int *)realloc(ipvl,memvl*2*sizeof(double));
    }
#ifdef ATOMLINEAR
    j = 1; nj = nmol;
    if(ipair==2){
      nj=nmol/2+1;
      if((nmol%2)==0 && (i>=nmol/2)) nj = nmol/2;
    }
#else  // good old i<j
    j = ((ipair==2)?i+1:0); 
    nj = nmol;
#endif
    for(kk=0;j<nj;j++){
#ifdef ATOMLINEAR
      int jj =  ((j+i)%nmol);
#else 
      int jj = j;
#endif
      if(!exclude(i,jj)) {
	double dis2=0.;
	for(k=0;k<DIM;k++) {double dx = xdis[k+kk];dis2 += dx*dx;}
	if(dis2<skin2){ipvl[2*nvl] = i;ipvl[2*nvl+1] = jj; nvl++;}
	kk += DIM;
      }
    }
    llist[i] = nvl;
  }
  delete[] xdis;
}


/*! 
  checks old positions agains new posistions to see if 
  neigbor list needs to be updated
*/
void NEIGHBOR :: checkSkin(COORDS &coords)
{
  static double *xo=NULL;
  int i,k,is,ie,ir=0;
  double rd,xdis[3];
  double *x=coords.x();

  double skin2 = skin*skin/4.;
  decomp1d(nmol,size,rank,&is,&ie);
  if(xo==NULL){
    xo = new double[(ie-is)*DIM];
    ir = 1; memvl = nmol;
    ipvl = (int *)malloc(memvl*2*sizeof(double));
  }
  for(i=is;i<ie && ir==0;i++){
    int ii = i*DIM;  
    int jj = (i-is)*DIM;
    for(k=0;k<DIM;k++) xdis[k] = x[ii+k] - xo[jj+k];
    for(rd=0.,k=0;k<DIM;k++) rd += xdis[k]*xdis[k];
    if(rd>=skin2) ir=1;
  }
#ifdef PARA
  {
    int ii;
    MPI_Allreduce(&ir,&ii,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    ir=ii;
  /*    MPI_Reduce(&ir,&ii,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD); */
  /* MPI_Bcast(&ir,1,MPI_INT,0,MPI_COMM_WORLD);  */
  }
#endif
  
  if(ir==1){
    ie = (ie-is)*DIM;
    for(i=0;i<ie;i++) xo[i] = x[i+is*DIM]; // store posistion
    getVlist(x);
  }
}

/*-------------------------------------------------------------*/
/*! get the next neighbors in the list using the no lists O(N^2) */
int NEIGHBOR :: nolist(void)
{
  static int inow,jnow,nj;
  int ni,i2,ibegin,iend;
  
  if(nintr==-2) {fprintf(stderr,"Why are we here nolist?\n");return 1;}
#ifdef ATOMLINEAR
  decomp1d(nmol,size,rank,&ibegin,&iend);
#else
  decomp_trig(nmol,size,rank,&ibegin,&iend);
#endif
  if(nintr==-1) {inow=ibegin;jnow=-1;}
  ni = 0;
  nnow = nix;
  nintr += nix;
  for(;inow<iend;inow++){
#ifdef ATOMLINEAR
    if(jnow==-1) jnow=0; // j>i 
    nj = nmol;
    if(ipair==0){
      nj=nmol/2+1;
      if((nmol%2)==0 && (inow>=nmol/2)) nj = nmol/2;
    }
#else
    if(jnow==-1 && ipair==0) jnow=inow; // j>i 
    nj = nmol;
#endif
    for(++jnow;jnow<nj;jnow++){
#ifdef ATOMLINEAR
      int jn = (jnow+inow)%nmol;
#else
      int jn = jnow;
#endif
      if(!exclude(jn,inow)){ // make sure i!=j
	i2 = ni*2;  ip[i2] = inow;  ip[i2+1] = jn;
	ni++;
	if(ni==nix) return 0; // not done
      }
    }
    jnow = -1;
  }
  
  nnow = ni;
  nintr = -2;
  return 1; // done
}
