
/*! \file ngbrplink.cc
    \brief file that contains function for neigbor lists using link lists and parallel processors
    
    \func servReq() sends coordinates out
    \func recvReq() recieve coordinates from processors
    \func 
    
*/

#include "mdmol.h++"

//--------------------------------------------------------

#ifdef PARA
void NEIGHBOR :: servReq(double x[]){
  const int NINTREQ=3;
  int flag;
  MPI_Status status;
  static MPI_Request req[2];
  do { // processes any return data
    MPI_Iprobe(MPI_ANY_SOURCE,TAGREQ,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      int ireq[NINTREQ];
      int proc = status.MPI_SOURCE;
      MPI_Recv(ireq,NINTREQ,MPI_INT,proc,TAGREQ,MPI_COMM_WORLD,&status);
      int jc = ireq[0];
      int nj = nlist[jc];      
      int *isbuf = new int[nj+NINTREQ];
      isbuf[nj] = ireq[0]; isbuf[nj+1] = ireq[1]; isbuf[nj+2] = ireq[2];
      if(nj>0){
	double *dbuf = new double[DIM*nj];
	getLink(jc,isbuf);
	getMyvals(nj,x,dbuf,isbuf);
	MPI_Isend(dbuf,nj*DIM,MPI_DOUBLE,proc,TAGSERVD,MPI_COMM_WORLD,&req[0]);
	MPI_Wait(&req[0],&status);
	delete[] dbuf;
      }
      MPI_Isend(isbuf,nj+NINTREQ,MPI_INT,proc,TAGSERVI,MPI_COMM_WORLD,&req[1]);
      MPI_Wait(&req[1],&status);
      delete[] isbuf;
    }
  } while (flag);
}

/*-----------------------------------------------------------------*/
int NEIGHBOR :: recvReq(int ipend[],double x[]){
  int flag,irec;
  MPI_Status status;
  irec = nintr = nnow = 0;
  do { // processes any return data
    MPI_Iprobe(MPI_ANY_SOURCE,TAGSERVI,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      irec = 1;
      int num,numt,ic,jc;
      int req = status.MPI_SOURCE;
      MPI_Get_count(&status,MPI_INT,&num);      
      int *jbuf = new int[num];
      MPI_Recv(jbuf,num,MPI_INT,req,TAGSERVI,MPI_COMM_WORLD,&status);
      ipend[jbuf[num-1]] += 1; // tag we received this message
      ic = jbuf[num-2];
      jc = jbuf[num-3];
      num -= 3; // remove icell,jcell and tag from buffer
      if(num != 0){
	MPI_Probe(req,TAGSERVD,MPI_COMM_WORLD,&status); // get particle info
	MPI_Get_count(&status,MPI_DOUBLE,&numt);    
	if(numt != num*DIM){ // need to tag individual message incase of async?
	  fprintf(stderr,"ERROR!!!! # don't match recieve cell part\n");
	  fprintf(stderr,"sent=%d != now=%d (rank=%d)\n",num*DIM*2,numt,rank);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
	double *rbuf = new double[numt];
	MPI_Recv(rbuf,numt,MPI_DOUBLE,req,TAGSERVD,MPI_COMM_WORLD,&status);
	int ni = nlist[ic];
	int nj = num;
	int *isbuf = new int[ni];
	double *xr = new double[DIM*ni];
	getLink(ic,isbuf);
	getMyvals(ni,x,xr,isbuf);
	if(nix < ni*nj){ 
	  fprintf(stderr,"ERROR: in link list, not enough space\n");
	  fprintf(stderr,"\tnnow=%d < ni(%d)*nj(%d)=%d\n",nix,ni,nj,ni*nj);
	  exit(1);
	}
	nnow = packxi(ni,nj,isbuf,jbuf,xr,rbuf,ip,xi);
	delete[] isbuf;
	delete[] xr;
	delete[] rbuf;
      }
      delete[] jbuf;
      if(nnow!=0) return 0;
    }
  } while (flag); // never gets here ;-) we could try to receive more mess.
  if(irec==1) return 0;
  return 1;
}
#endif

//--------------------------------------------------------
#ifdef PARA
void NEIGHBOR :: sendcoords(int np,int ipm[],int jcell,double x[],double v[]){
  static MPI_Request rq[2];
  MPI_Status status;
  int i,k,j,ip;
  double *sndbuf = new double[np*DIM*2];

  ip = ncell*ncell*ncell;
  int proc = icpmap[jcell]; // proc to send part to
  for(i=0;i<np-1;i++){ // pack up coords
    j = i*DIM*2; ip = ipm[i];
    for(k=0;k<DIM;k++) sndbuf[j+k] = x[ip+k];
    for(k=0;k<DIM;k++) sndbuf[j+k + DIM] = v[ip+k];
  }
  // send long message first
  MPI_Isend(sndbuf,(np-1)*DIM*2,MPI_DOUBLE,proc,TAGUPDD,MPI_COMM_WORLD,&rq[0]);
  MPI_Isend(ipm,np,MPI_INT,proc,TAGUPDI,MPI_COMM_WORLD,&rq[1]);  
  MPI_Wait(&rq[0],&status);
  MPI_Wait(&rq[1],&status);
  delete[] sndbuf;
}

//--------------------------------------------------------
void NEIGHBOR :: recvcoords(double x[],double v[],int proc[]){
  int flag;
  MPI_Status status;
  do { // processes any return data
    MPI_Iprobe(MPI_ANY_SOURCE,TAGUPDI,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      int num,numt,i;
      int req = status.MPI_SOURCE;
      MPI_Get_count(&status,MPI_INT,&num);      
      int *ipmapb = new int[num];
      MPI_Recv(ipmapb,num,MPI_INT,req,TAGUPDI,MPI_COMM_WORLD,&status);
      proc[ipmapb[num-1]] = 1; // tag we received this message
      num--; // remove tag from count
      for(i=0;i<num;i++) ipmap[nmp+i] = ipmapb[i];
      nmp += num; // store received particles number
      MPI_Probe(req,TAGUPDD,MPI_COMM_WORLD,&status); // get particle info
      MPI_Get_count(&status,MPI_DOUBLE,&numt);     
      if(numt != num*DIM*2){ // need to tag individual message incase of async?
	fprintf(stderr,"ERROR!!!! # don't match recieve coords\n");
	fprintf(stderr,"sent=%d != now=%d\n",num*DIM*2,numt);
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      double *rbuf = new double[numt];
      MPI_Recv(rbuf,numt,MPI_DOUBLE,req,TAGUPDD,MPI_COMM_WORLD,&status);
      // unpack 
      for(i=0;i<num;i++){
	int j,k,ip; j = i*DIM*2; ip = ipmapb[i];
	for(k=0;k<DIM;k++) x[ip+k] = rbuf[j+k];
	for(k=0;k<DIM;k++) v[ip+k] = rbuf[j+k + DIM];
      }
      delete[] ipmapb;
      delete[] rbuf;
    }
  } while (flag);
}

/* -------------- collect particles to master node --------------------*/

void NEIGHBOR :: recvcdsm(double x[],double v[],int proc[]){
  int flag;
  MPI_Status status;
  do { // processes any return data
    MPI_Iprobe(MPI_ANY_SOURCE,TAGCMI,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      int num,numt;
      int req = status.MPI_SOURCE;
      MPI_Get_count(&status,MPI_INT,&num);      
      int *ipmapb = new int[num];
      MPI_Recv(ipmapb,num,MPI_INT,req,TAGCMI,MPI_COMM_WORLD,&status);
      proc[req] += 1; // tag we received this message
      
      MPI_Probe(req,TAGCMD,MPI_COMM_WORLD,&status); // get particle info
      MPI_Get_count(&status,MPI_DOUBLE,&numt);     
      if(numt != num*DIM*2){ // need to tag individual message incase of async?
	fprintf(stderr,"ERROR!!!! # don't match collect part\n");
	fprintf(stderr,"sent=%d != now=%d\n",num*DIM*2,numt);
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      double *rec = new double[numt];
      MPI_Recv(rec,numt,MPI_DOUBLE,req,TAGCMD,MPI_COMM_WORLD,&status);
      scatMyvals(num,x,rec,ipmapb);       // unpack 
      scatMyvals(num,v,rec+num*DIM,ipmapb);       // unpack 
      delete[] ipmapb;
      delete[] rec;
    }
  } while (flag);
}
#endif

