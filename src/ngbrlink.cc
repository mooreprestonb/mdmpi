
/*! \file ngbrlink.cc
    \brief file that contains function for neigbor lists using link lists
    
    
*/

#include "mdmol.h++"

/*--------------------------------------------------------*/
int getCngbr(int jcell,int icell,int ncell,int iperd)
{
  int k,ncell2,ic[3],jc[3];

  ncell2 = ncell*ncell;
  
  ic[0] = icell/ncell2;               /* get cell i,j,k index */
  ic[1] = (icell - ic[0]*ncell2)/ncell;
  ic[2] = (icell - ic[0]*ncell2 - ic[1]*ncell);

  switch(jcell){
  case 0:  jc[0]=ic[0]  ; jc[1]=ic[1]  ; jc[2]=ic[2]  ; break;
  case 1:  jc[0]=ic[0]+1; jc[1]=ic[1]  ; jc[2]=ic[2]  ; break;
  case 2:  jc[0]=ic[0]-1; jc[1]=ic[1]+1; jc[2]=ic[2]  ; break;
  case 3:  jc[0]=ic[0]  ; jc[1]=ic[1]+1; jc[2]=ic[2]  ; break;
  case 4:  jc[0]=ic[0]+1; jc[1]=ic[1]+1; jc[2]=ic[2]  ; break;
  case 5:  jc[0]=ic[0]-1; jc[1]=ic[1]-1; jc[2]=ic[2]+1; break;
  case 6:  jc[0]=ic[0]  ; jc[1]=ic[1]-1; jc[2]=ic[2]+1; break;
  case 7:  jc[0]=ic[0]+1; jc[1]=ic[1]-1; jc[2]=ic[2]+1; break;
  case 8:  jc[0]=ic[0]-1; jc[1]=ic[1]  ; jc[2]=ic[2]+1; break;
  case 9:  jc[0]=ic[0]  ; jc[1]=ic[1]  ; jc[2]=ic[2]+1; break;
  case 10: jc[0]=ic[0]+1; jc[1]=ic[1]  ; jc[2]=ic[2]+1; break;
  case 11: jc[0]=ic[0]-1; jc[1]=ic[1]+1; jc[2]=ic[2]+1; break;
  case 12: jc[0]=ic[0]  ; jc[1]=ic[1]+1; jc[2]=ic[2]+1; break;
  case 13: jc[0]=ic[0]+1; jc[1]=ic[1]+1; jc[2]=ic[2]+1; break;
    // bottom half if you arn't using newton's 3rd law    
  case 14: jc[0]=ic[0]-1; jc[1]=ic[1]  ; jc[2]=ic[2]  ; break;
  case 15: jc[0]=ic[0]+1; jc[1]=ic[1]-1; jc[2]=ic[2]  ; break;
  case 16: jc[0]=ic[0]  ; jc[1]=ic[1]-1; jc[2]=ic[2]  ; break;
  case 17: jc[0]=ic[0]-1; jc[1]=ic[1]-1; jc[2]=ic[2]  ; break;
  case 18: jc[0]=ic[0]+1; jc[1]=ic[1]+1; jc[2]=ic[2]-1; break;
  case 19: jc[0]=ic[0]  ; jc[1]=ic[1]+1; jc[2]=ic[2]-1; break;
  case 20: jc[0]=ic[0]-1; jc[1]=ic[1]+1; jc[2]=ic[2]-1; break;
  case 21: jc[0]=ic[0]+1; jc[1]=ic[1]  ; jc[2]=ic[2]-1; break;
  case 22: jc[0]=ic[0]  ; jc[1]=ic[1]  ; jc[2]=ic[2]-1; break;
  case 23: jc[0]=ic[0]-1; jc[1]=ic[1]  ; jc[2]=ic[2]-1; break;
  case 24: jc[0]=ic[0]+1; jc[1]=ic[1]-1; jc[2]=ic[2]-1; break;
  case 25: jc[0]=ic[0]  ; jc[1]=ic[1]-1; jc[2]=ic[2]-1; break;
  case 26: jc[0]=ic[0]-1; jc[1]=ic[1]-1; jc[2]=ic[2]-1; break;
  default: fprintf(stderr,"Internal ERROR in getCngbr!\n");exit(1); break;
  }
  
  /* apply boundry conditions */  
  for(k=0;k<3;k++) if(jc[k]==ncell) if(iperd>k) jc[k]=0; else return -1;
  for(k=0;k<3;k++) if(jc[k]== -1) if(iperd>k) jc[k]=ncell-1; else return -1;
  
  /* get jcell from i,j,k */
  for(jcell=jc[0],k=1;k<3;k++) jcell = jcell*ncell+jc[k];/*ic*n2+jc*n+kc*/
  return jcell;
}

/*-----------------------------------------------------------------*/
void NEIGHBOR :: getLnklist(double x[],int npnow)
{
  int i,k,ic[3],icell;
  
  int ncellm1 = ncell-1;
  for(i=0;i<nmol;i++) llist[i] = -1;
  k = ncell*ncell*ncell;
  for(i=0;i<k;i++) hlist[i] = -1;
  for(i=0;i<k;i++) nlist[i] = 0;

  double *sx = new double[npnow*DIM];
  for(i=0;i<npnow*DIM;i++) sx[i] = x[i];
  iperiodic(npnow,sx);   // reduce to 1/cell
  for(i=npnow-1;i>=0;i--){
    int ii = i*DIM;
    for(k=0;k<DIM;k++) ic[k] = (int)((sx[k+ii]+.5)*ncell);
    for(k=0;k<DIM;k++) ic[k] = MAX(ic[k],0);  // make sure we are not neg
    for(k=0;k<DIM;k++) ic[k] = MIN(ic[k],ncellm1); // make sure we arn't over
    for(icell=ic[0],k=1;k<DIM;k++) icell = icell*ncell+ic[k];/*ic*n2+jc*n+kc*/
    llist[i] = hlist[icell];
    hlist[icell] = i;
  }
  delete[] sx;

  for(i=0;i<ncell*ncell*ncell;i++) {
    int n=0;
    int ii = hlist[i];
    while(ii != -1){ n++;ii = llist[ii];}
    nlist[i] = n;
  }
#ifdef DEBUG
  printf("Link head\n");
  for(i=0;i<ncell*ncell*ncell;i++) 
    printf("%d/hist %d %d\n",rank,i,hlist[i]);
  printf("Link List\n");
  for(i=0;i<npnow;i++) 
    printf("%d/llist %d %d %d\n",rank,i,llist[i],ipmap[i]);
#endif
}

/*-----------------------------------------------------------------*/
int NEIGHBOR :: getLink(int ic,int ip[])
{
  int nn=0;int jj = hlist[ic];
  while(jj != -1){ip[nn++] = jj;jj = llist[jj];} // end while(jj!=-1)
  return nn;
}

/*-----------------------------------------------------------------*/
int NEIGHBOR :: packxi(int ni,int nj,int ibuf[],int jbuf[],
		       double xib[] ,double xjb[],int ip[],double xi[])
{
  int i,j,k;
  int nnow=0;
  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      if(!exclude(ibuf[i],jbuf[j])){
	int id = ibuf[i]; int jd = jbuf[j];
	int nn = 2*nnow; ip[nn] = id; ip[nn+1] = jd;
	nn *= DIM;
	for(k=0;k<DIM;k++){xi[nn+k]=xib[i*DIM+k];xi[nn+k+DIM]=xjb[j*DIM+k];}
	nnow++;
      }
    }
  }
  return nnow;
}
/*-----------------------------------------------------------------*/
int NEIGHBOR :: packxi(int ni,int nj,int ibuf[],int jbuf[],int ip[])
{
  int i,j;
  int nnow=0;
  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      if(!exclude(ibuf[i],jbuf[j])){
	int id = ibuf[i]; int jd = jbuf[j];
	int nn = 2*nnow; ip[nn] = id; ip[nn+1] = jd;
	nnow++;
      }
    }
  }
  return nnow;
}
/*-----------------------------------------------------------------*/
int NEIGHBOR :: interLink(void)
{
  static int *ipend=NULL,*nproc;
  static int ic=0,js=0,ii,done,sdone;
  int i,nn,nil,nih,jc;
#ifdef PARA 
  static MPI_Request *req;
#endif

  int nnp = (ipair==4?14:27);
  if(ipend==NULL) {
    ipend = new int[nnp];
    nproc = new int[size];
#ifdef PARA
    req = new MPI_Request[size];
#endif
  }

  if(nintr==-2) {fprintf(stderr,"why are we here interLink\n");return 1;}
  nn = ncell*ncell*ncell;
  decomp1d(nn,size,rank,&nil,&nih); 
  if(nintr == -1) { // initialize static variables 
    ic = nil; js = 0;ii = -1;done=0;sdone=0;
    for(i=0;i<size;i++) nproc[i] = -1;
    for(i=0;i<nnp;i++) ipend[i] = 0; // zero all cells
  }

  nintr = nnow = 0;
  for(;ic<nih;ic++){ // loop until we are done with our cell 
    if(nnow != 0) return 0;
#ifdef PARAO
    servReq(x);
    if(recvReq(ipend,x)==0) return 0;
#endif
    int ni = nlist[ic];
    if(ni==0) {
      ic++;
      for(i=0;i<nnp;i++) ipend[i] += 1;
      return 0;
    }
    int *isbuf = new int[ni];
    getLink(ic,isbuf);
    // double *xr = new double[DIM*ni];
    // getMyvals(ni,x,xr,isbuf);
    
    // same cell
    if(js==0){
      if(nix < ni*ni){ // compute stuff?
	fprintf(stderr,"ERROR: in link list, not enough space\n");
	fprintf(stderr,"\tnnow=%d < ni(%d)^2=%d\n",nix,ni,ni*ni);
	exit(1);
      }
      for(i=0;i<ni;i++){
	int j = ((ipair==4)?i+1:0);
	for(;j<ni;j++){
	  if(!exclude(i,j)){
	    int id = isbuf[i]; int jd = isbuf[j];
	    int nn = 2*nnow; ip[nn] = id; ip[nn+1] = jd;
	    nnow++;
	  }
	}
      }
      ipend[js] += 1;
      js++;
      delete[] isbuf;
      // delete[] xr;
      return 0;
    }
    // different cells 
    for(;js<nnp;js++){
#ifdef PARAO
      servReq(x);
      if(recvReq(ipend,x)==0) {
	delete[] isbuf;
	delete[] xr;
	return 0;
      }
#endif
      jc = getCngbr(js,ic,ncell,iperd);
      if(jc != -1){
	if(icpmap[jc] == rank){
	  int nj = nlist[jc];
	  ipend[js] += 1;
	  if(nj != 0){
	    int *jsbuf = new int[nj];
	    // double *xj = new double[DIM*nj];
	    getLink(jc,jsbuf);
	    // getMyvals(nj,x,xj,jsbuf);	  
	    if(nix < ni*nj){ // compute stuff?
	      fprintf(stderr,"ERROR: in link list, not enough space\n");
	      fprintf(stderr,"\tnnow=%d < ni(%d)*nj(%d)=%d\n",nix,ni,nj,ni*nj);
	      exit(1);
	    }
	    nnow = packxi(ni,nj,isbuf,jsbuf,ip);
	    js++;
	    delete[] jsbuf;
	    delete[] isbuf;
	    return 0;
	  }
	} else {  // we need to request these particles
#ifdef PARA
	  int ireq[3];
	  MPI_Status status;
	  MPI_Request req;
	  ireq[0] = jc; ireq[1] = ic; ireq[2] = js;
	  MPI_Isend(ireq,3,MPI_INT,icpmap[jc],TAGREQ,MPI_COMM_WORLD,&req);
	  MPI_Wait(&req,&status);
#endif
	}
      } else ipend[js] += 1; // end if jc != -1
    } // for js
    js = 0;
    delete[] isbuf;
    // delete[] xr;
  } 
#ifdef PARAO
  int itr=0;
  while(!done){
    servReq(x);
    doneReq(nproc);
    if(alldone(nnp,ipend,nih-nil)){ // are we done? 
      if(!sdone){ // is everybody done?
	for(int p=0;p<size;p++){
	  if(p==rank) nproc[p] = p;
	  else MPI_Isend(NULL,0,MPI_INT,p,DONETAG,MPI_COMM_WORLD,&req[p]);
	} 
	sdone = 1;
      }
    } else {
      if(recvReq(ipend,x)==0) return 0;
    }
    if(itr>ITERMAX){
      fprintf(stderr,"WARNING: ITERMAX REACHED INTRL rank=%d - %d, %d,%d\n",
	      rank,itr,alldone(nnp,ipend,nih-nil),nih-nil);
      for(i=0;i<nnp;i++) fprintf(stderr,"proc[%d] = %d\n",i,ipend[i]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    for(i=0;i<size;i++) if(nproc[i]!=i) break;
    if(i==size) done=1;
    itr++;
  }
  if(noReq(rank)) {
    fprintf(stderr,"Some requests not done in INTRL?\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
#endif
  nintr = -2;
  return 1;
}
/* -------------- move any particles that we need to --------------------*/
void rmipmap(int nmp,int ipmap[],int ii,int rank){
  int i;
  for(i=0;i<nmp;i++){
    if(ipmap[i] == ii){ipmap[i] = -1;break;}
  }
  if(i==nmp){
    fprintf(stderr,"ERROR: in rmipmap, ii=%d, rank=%d\n",ii,rank);
#ifdef PARA
    MPI_Abort(MPI_COMM_WORLD,1);    
#endif
    exit(1);
  }
}

/* -------------- move any particles that we need to --------------------*/
void NEIGHBOR :: updateLink(COORDS &coords){
#ifdef PARA
  if(noReq(rank)) {
    fprintf(stderr,"Some requests not done in entering Update?\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int i,ic,js,np,start,end;
  scatMyvals(nmp,coords.x(),xmpart,ipmap);
  scatMyvals(nmp,coords.v(),xmpart+nmp*DIM,ipmap);
  getLnklist(xmpart,nmp);
  np = ncell*ncell*ncell;
  decomp1d(np,size,rank,&start,&end); 
  int *ipmapb = new int[nmol];
  int nnp = (ipair==4?14:27);
  int *ipend = new int[nnp];
  int nc3 = (end-start);
  for(i=0;i<nnp;i++) ipend[i] = 0;
  for(np=0,ic=start;ic<end;ic++){
    for(js=0;js<nnp;js++){
      recvcoords(coords.x(),coords.v(),ipend); //recieve any
      int jc = getCngbr(js,ic,ncell,iperd);
      if((jc != -1) && (jc<start || jc >= end)){ // not our cell so we send
	int ii = hlist[jc];
	while(ii != -1){
	  rmipmap(nmp,ipmap,ii,rank); // remove from ipmap bit
	  ipmapb[np++] = ii; 
	  ii = llist[ii];
	}
	if(np != 0){ // we have np particles to send
	  ipmapb[np++] = js;
	  sendcoords(np,ipmapb,js,coords.x(),coords.v());
	} else ipend[js] += 1; // mark that we have sent/recieved it
      } else {
	ipend[js] += 1; // mark that we have recieved it
      }
    }
  }
  // receive part 
  int itr = 0;
  while(!alldone(nnp,ipend,nc3)){ // are we done? 
    recvcoords(coords.x(),coords.v(),ipend); //get rest mess
    if(itr>ITERMAX){
      fprintf(stderr,"WARNING: ITERMAX REACHED UPD rank=%d - %d, %d,%d\n",
	      rank,itr,alldone(nnp,ipend,nc3),nc3);
      for(i=0;i<size;i++) fprintf(stderr,"proc[%d] = %d\n",i,ipend[i]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    itr++;
  }

  // loop over part and get ours!
  int nmpt=0;
  for(i=0;i<nmp;i++) if(ipmap[i] != -1) ipmap[nmpt++] = ipmap[i];
  nmp = nmpt;
  if(nmp >= nmpm) {
    nmpm = (nmp+MAXALLOC);  
    xmpart= (double *)realloc(xmpart,nmpm*DIM*3*sizeof(double));  
  }
  getMyvals(nmp,coords.x(),xmpart,ipmap);
  getMyvals(nmp,coords.v(),xmpart+nmp*DIM,ipmap);
  
  delete[] ipmapb;
  delete[] ipend;

  if(noReq(rank)) {
    fprintf(stderr,"Some requests not done in leaving Update?\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
#endif
  getLnklist(xmpart,nmp);
}
//---------------------------------------------------------------------
