
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

/*-------------------------------------------------------------*/

void decomp1d(int n,int size,int rank,int *s,int *e)
{
  int nlocal, deficit;

  nlocal= n/size;     /*nice local size */
  *s=rank * nlocal;
  deficit=n%size;     /* remainder */
  *s += ((rank < deficit) ? rank : deficit); /* add deficit to s if nec*/
  if (rank < deficit) nlocal++;
  *e = *s + nlocal;            /* calculate end point */
  if (*e > n || rank == size-1) *e = n;
}

#define WARNTRIG
/*-------------------------------------------------------------*/
void decomp_trigi_old(int n,int size,int rank,int *ibegin,int *iend)
{ 
  unsigned int i;
  int s,e;

  if(n>((INT_MAX/n)/4)){
    /* npairs to big for integer so we spread i evenly */
#ifdef WARNTRIG
    fprintf(stderr,"Warning: npairs (%d) is too big for integer\n",n);
    fprintf(stderr,"\t try using different neighbor lists\n");
    fprintf(stderr,"\t %d %d\n",n,INT_MAX);
#endif
    decomp1d(n,size,rank,ibegin,iend);
  } else {
    i = n*(n-1)/2;
    decomp1d(i,size,rank,&s,&e);
    i = 1+4*n*(n-1);
    *ibegin = (int)((-1.+2.*n-sqrt((double)(i-8*s)))/2.);
    *iend = (int)((-1.+2.*n-sqrt((double)(i-8*e)))/2.);
  }
}
/*-------------------------------------------------------------*/
void decomp_trig(int n,int size,int rank,int *ibegin,int *iend)
{ 
  unsigned int nlocal,nl;
  int s,e;

  if(n>(INT_MAX/n)){
    /* npairs to big for integer so we spread i evenly */
#ifdef WARNTRIG
    fprintf(stderr,"Warning: npairs (%d) is too big for integer 2\n",n);
    fprintf(stderr,"\t try using different neighbor lists\n");
    fprintf(stderr,"\t %d %d\n",n,INT_MAX);
#endif
    decomp1d(n,size,rank,ibegin,iend);
  } else {
    nl = (n%2) ? n*((n-1)/2) : (n/2)*(n-1);
    nlocal = nl/size;
    rank = size-1-rank;
    s = rank * nlocal;
    e = s + nlocal;
    if(rank==size-1) *ibegin=0;
    else *ibegin = n-(int)(sqrt((double)(2*e+.25))+.5);
    if(rank==0) *iend = n;
    else  *iend = n-(int)(sqrt((double)(2*s+.25))+.5);
  }
}
/*-------------------------------------------------------------*/

