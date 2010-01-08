/*! \file ngbrpara.cc 
  this file containes neighbor list functions for dealing
  with paralization of the neighbor lists and particles across
  processors
*/
#include "mdmol.h++"

#ifdef PARA
#include <mpi.h>
#endif

//! all the process have received all the need to
int alldone(int nr,int iproc[],int n){
  int done=1; 
  for(int i=0;i<nr;i++) if(iproc[i] != n){done=0;break;}
  return done;
}
//! make sure that there are not outstanding messaged to be processed
int noReq(int rank)
{
  int ierr=0;
#ifdef PARA
  int flag;
  MPI_Status status;
  MPI_Barrier(MPI_COMM_WORLD);  // we should all be at the same place!
  do { // serve any requests
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      int num;
      int req = status.MPI_SOURCE;
      MPI_Get_count(&status,MPI_CHAR,&num);      
      char *cbuf = new char[num];
      MPI_Recv(cbuf,num,MPI_CHAR,req,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      ierr++;
      fprintf(stderr,"Messages %d not rec rank=%d from %d! num=%d, tag=%d\n",
	      ierr,rank,status.MPI_SOURCE,num,status.MPI_TAG);
      delete[] cbuf;
    }
  } while (flag);
  if(ierr){
    fprintf(stderr,"There were %d messages not proccessed\n",ierr);
    fflush(stderr);  fflush(stdout); sleep(2);  
  }
#endif
  return ierr;
}
//! done request... all the request are done
void doneReq(int nproc[])
{
#ifdef PARA
  int flag;
  MPI_Status status;
  do { // serve any requests
    MPI_Iprobe(MPI_ANY_SOURCE,DONETAG,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      int req = status.MPI_SOURCE;
      MPI_Recv(NULL,0,MPI_INT,req,DONETAG,MPI_COMM_WORLD,&status);
      nproc[req] = req;
    }
  } while (flag);
#endif
}
//! scaters x into xp with map ipmap
void scatMyvals(int npart,double x[],double xp[],int ipmap[])
{
  int i,k,is,ii;
  
  for(i=0;i<npart;i++) {
    ii = i*DIM;    is=ipmap[i]*DIM;
    for(k=0;k<DIM;k++) x[is+k] = xp[ii+k];
  }
}
//! gathers x from xp with map ipmap 
void getMyvals(int npart,double x[],double xp[],int ipmap[])
{
  int i,k,is,ii;
  
  for(i=0;i<npart;i++) {
    ii = i*DIM;    is=ipmap[i]*DIM;
    for(k=0;k<DIM;k++) xp[ii+k] = x[is+k];
  }
}
/*! initializes myparticales  */
double * NEIGHBOR :: initMypart(COORDS &coords)
{
  int i,start,end,np;
  np = 0;
  counts = new int[size];
  disp = new int[size];
  
  switch(ipair){ /* initial neigbor lists */
  case 0:  case 1: case 2: case 3: // no and verlet lists continuous npart
    decomp1d(nmol,size,rank,&start,&end);
    np = end-start;
    ipmap = new int[np];
    for(i=0;i<np;i++) ipmap[i] = i+start; // map total part to mypart
    break;
  case 4: case 5:
    np = ncell*ncell*ncell;
    decomp1d(np,size,rank,&start,&end); 
    ipmap = new int[nmol];
    getLnklist(coords.x(),nmol);
    for(np=0,i=start;i<end;i++){
      int ii = hlist[i];
      while(ii != -1){
	ipmap[np++] = ii; 
	ii = llist[ii];
      }
    }
    break;
  default:
    fprintf(stderr,"ERROR: in InitMyPart ipair=%d\n",ipair); exit(1);
  }
  
  for(i=0;i<size;i++) disp[i]=0;
  disp[rank]=nmp=np;
#ifdef PARA
  nmpm = (np+1);
  MPI_Allreduce(disp,counts,size,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  xmpart = (double *)realloc(xmpart,nmpm*DIM*3*sizeof(double));  
#else
  xmpart = coords.x();
  nmpm = (np+1);
  counts[rank]=disp[rank];
#endif  
  counts[0] *= DIM;  disp[0] = 0;
  for(i=1;i<size;i++){
    counts[i] *= DIM;
    disp[i] = disp[i-1]+counts[i-1];
  }

  return xmpart;
}

/*! get the data onto the centrial processor */
void NEIGHBOR :: getcm(double x[],double v[]){
  switch(ipair){
  case 0: case 1: case 2: case 3: break; // replicated data x is already here
  case 4: case 5: 
#ifdef PARA
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
      int i;
      int *ipend = new int[size];
      for(i=0;i<size;i++) ipend[i]=0;
      scatMyvals(nmp,x,xmpart,ipmap);
      scatMyvals(nmp,v,xmpart+nmp*DIM,ipmap);
      ipend[0] += 1;
      int itr=0;
      while(!alldone(size,ipend,1)){ // are we done? 
	recvcdsm(x,v,ipend); //get more messages
	if(itr>ITERMAX){
	  fprintf(stderr,"WARNING: ITERMAX REACHED getcm rank=%d-%d, %d,%d\n",
		  rank,itr,alldone(size,ipend,1),1);
	  for(i=0;i<size;i++) fprintf(stderr,"proc[%d] = %d\n",i,ipend[i]);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
	itr++;
      }
      delete[] ipend;
    } else {
      MPI_Status status;
      MPI_Request req[2];
      MPI_Isend(ipmap,nmp,MPI_INT,0,TAGCMI,MPI_COMM_WORLD,&req[1]);  
      MPI_Isend(xmpart,nmp*DIM*2,MPI_DOUBLE,0,TAGCMD,MPI_COMM_WORLD,&req[0]);
      MPI_Wait(&req[1],&status);
      MPI_Wait(&req[0],&status);
    }
    if(noReq(rank)) {
      fprintf(stderr,"Some requests not done in getcm\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
#endif
    break;
  default:fprintf(stderr,"ERROR: in get_neigbors ipair=%d\n",ipair); exit(1);
  }
}
/*--------------------------------------------------------*/
