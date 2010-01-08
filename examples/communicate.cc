
// test iprobe and recieve and etc....

#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

enum TAGS {SENDTAG,RECVTAGI,RECVTAGD,DONETAG};

#define MAX(A,B) (((A)>(B))?(A):(B))

//--------------------------------------------------------
void noReq(int rank)
{
  int flag,ierr=0;
  MPI_Status status;
  do { // serve any requests
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      ierr++;
      fprintf(stderr,"Messages from %d to %d not done!\n",
	      status.MPI_SOURCE,rank);
    }
  } while (flag);
  if(ierr!=0){
    fprintf(stderr,"There were %d messages not proccessed\n",ierr);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}
//--------------------------------------------------------
void donReq(int nproc[])
{
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
}
//--------------------------------------------------------
void serReq(int rank,int size)
{
  int flag;
  static MPI_Request request[2]; static MPI_Status status;

  do { // serve any requests
    MPI_Iprobe(MPI_ANY_SOURCE,SENDTAG,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      int j, num, req = status.MPI_SOURCE;
      MPI_Get_count(&status,MPI_INT,&num);      
      if(num!=2) {
	printf("somebody (%d) requested %d numbers from %d?\n",req,num,rank);
	MPI_Abort(MPI_COMM_WORLD,1);
      } else {
	int rbuf[2];
	MPI_Recv(rbuf,2,MPI_INT,req,SENDTAG,MPI_COMM_WORLD,&status);
	// pakage
	double *array = (double *)malloc(rbuf[0]*sizeof(double)); //data
	// send
	MPI_Isend(array,rbuf[0],MPI_DOUBLE,req,RECVTAGD,MPI_COMM_WORLD,
		  &request[0]); // send long message first
	MPI_Isend(rbuf,2,MPI_INT,req,RECVTAGI,MPI_COMM_WORLD,&request[1]);
	free(array);
      }
    }
  } while (flag);
}
//--------------------------------------------------------
void resReq(int arproc[],int arnum[],int drproc[])
{
  int flag;
  MPI_Status status;
  do { // processes any return data
    MPI_Iprobe(MPI_ANY_SOURCE,RECVTAGI,MPI_COMM_WORLD,&flag,&status);
    if(flag){
      int num;
      int req = status.MPI_SOURCE;
      int rbuf[2];
      MPI_Get_count(&status,MPI_INT,&num);      
      if(num!=2){
	fprintf(stderr,"ERROR num!=2\n");MPI_Abort(MPI_COMM_WORLD,1);
      }
      MPI_Recv(rbuf,2,MPI_INT,req,RECVTAGI,MPI_COMM_WORLD,&status);
      num = rbuf[1];
      if(arproc[num] != req){
	fprintf(stderr,"ERROR!!!! no matching call %d %d != %d\n",
		num,arproc[num],req);
	MPI_Abort(MPI_COMM_WORLD,1);
      } else if (drproc[num] == 1){
	fprintf(stderr,"ERROR!!!! recieved call twice\n");
	MPI_Abort(MPI_COMM_WORLD,1);
      } else {
	drproc[num] = 1;
      }
      if(rbuf[0] != arnum[rbuf[1]]){ // 
	fprintf(stderr,"ERROR!!!! # of random numbers don't match #1\n");
	fprintf(stderr,"rbug[0]=%d != arnum[%d]=%d\n",
		rbuf[0],rbuf[1],arnum[rbuf[1]]);
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      // now that we know what is comming the get it with a blocking Recv.
      MPI_Probe(req,RECVTAGD,MPI_COMM_WORLD,&status);
      MPI_Get_count(&status,MPI_DOUBLE,&num);      
      if(rbuf[0] != num){
	fprintf(stderr,"ERROR!!!! # of random numbers don't match #2\n");
	fprintf(stderr,"rbug[0]=%d != num=%d\n",rbuf[0],num);
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      double *array = (double *)malloc(num*sizeof(double));
      MPI_Recv(array,rbuf[0],MPI_DOUBLE,req,RECVTAGD,MPI_COMM_WORLD,&status);
      average(rbuf[0],array); 
      free(array);
    }
  } while (flag);
}

//--------------------------------------------------------
int alldone(int nr,int drproc[]){
  int done=1; 
  for(int i=0;i<nr;i++) if(drproc[i] != 1){done=0;break;}
  return done;
}

//--------------------------------------------------------
void request(int rank,int size)
{
  static int *nproc = 0;
  static int *arproc,*drproc,*arnum;
  static MPI_Request *req;
  int i,n;
  if(nproc == 0) {
    n = MAXNUMS;
    nproc  = (int *)malloc(size*sizeof(int));
    arproc = (int *)malloc(3*n*sizeof(int));
    drproc = arproc + n;
    arnum  = arproc + 2*n;
    req = (MPI_Request *)malloc(size*sizeof(MPI_Request));
  }
  n = int(randme()*MAXNUMS);
  if(istep < 2) n = size;  
  switch(istep){
  case 0: // first pass all have the same numbers
    for(i=0;i<n;i++) {arproc[i] = i; arnum[i] = int(size);} break;
  case 1: // second we send different amounts...
    for(i=0;i<n;i++) {arproc[i] = i;arnum[i] = i+1;} break;
  default: // now we randomize
    for(i=0;i<n;i++) arproc[i] = int(randme()*size);
    for(i=0;i<n;i++) arnum[i] = int(randme()*MAXNUMS);
    break;
  }
  for(i=0;i<n;i++) drproc[i] = 0;
  for(i=0;i<size;i++) nproc[i] = -1;

  //for(i=0;i<n;i++)
  //    fprintf(fout,"rank=%d, %d %d %d\n",rank,i,arproc[i],arnum[i]);
  
  int done = 0, sdone = 0, itr = 0;
  while(!done){
    serReq(rank,size); // service requests
    if(itr<n) {  // we send request 
      if(arnum[itr]>0){
	int arbuf[2];
	arbuf[0] = arnum[itr];
	arbuf[1] = itr;
	MPI_Isend(arbuf,2,MPI_INT,arproc[itr],SENDTAG,MPI_COMM_WORLD,
		  &req[arproc[itr]]);
      } else {
	drproc[itr] = 1; // no need for requests :-)
      }
    }
    resReq(arproc,arnum,drproc); // receive requests 
    if(alldone(n,drproc)){ // are we done?
      int p;
      if(!sdone){ // is everybody done?
        for(p=0;p<size;p++){
          if(p==rank) nproc[p] = rank;
          else MPI_Send(NULL,0,MPI_INT,p,DONETAG,MPI_COMM_WORLD);
        }
        sdone = 1;
      } else {
        donReq(nproc);
        for(p=0;p<size;p++) if(nproc[p] != p) break;
        if(p==size) done = 1;
      }
    }
    if(itr>ITERMAX){
      fprintf(stderr,"WARNING: ITERMAX REACHED rank=%d - %d, %d %d %d\n",
	      rank,itr,sdone,alldone(n,drproc),done);
      for(i=0;i<size;i++) 
        fprintf(stderr,"proc[%d] = %d\n",i,nproc[i]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    ++itr;
  }

  if(rank==0) {
    fprintf(fout,"rank=%d, avg(%d)=%g, step=%d, itr=%d\n",
	    rank,n,average(0,(double *)nproc),istep,itr);
  }
  noReq(rank); // check to make sure we have no messages waiting
}
//--------------------------------------------------------
