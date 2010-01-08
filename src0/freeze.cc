// routines to freeze the certain atoms

#include <stdlib.h>
#include "freeze.h++"

const int DIM=3;
int nfreeze,*frlist=NULL;

void freeFreezeScr(void)
{
  if(frlist != NULL) free(frlist);
  nfreeze=0;
}

void setupFreeze(int n)
{
  int i;
  nfreeze=n;
  if(n>0){
    frlist = (int *)malloc(n*sizeof(int));
    for(i=0;i<n;i++) frlist[i] = i;
  }
}


void zeroFreeze(double a[])
{
  int i,k;
  for(i=0;i<nfreeze;i++){
    int ii=frlist[i]*DIM;
    for(k=0;k<DIM;k++) a[ii+k] = 0.;
  }
}
