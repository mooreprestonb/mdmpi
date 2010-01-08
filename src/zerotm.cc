/*! 
  \file zerotm.cc
  \brief Functions to zero the total momentum ie linear & angular 
  
  for now we only have the linear momentum and don't freeze atoms.
*/

/*! 
  \fn void zerotm(int n,int DIM,double v[],double m[])
  \brief Function to zero the total momentum 
  zeros the total momentum of velocity \a v with mass \a m

  \param n the number of elements 
  \param DIM the number of dimensions
  \param nfreeze the number of atoms frozen
  \param v vector of velocities with length \a DIM * \a n
  \param m vector of masses with length \a n
*/

#include "zerotm.h++"

void zerotm(int n,int DIM,int nfreeze,double v[],double *m)
{
  int i,k; 
  double tv[DIM]; 
  double rn = 0; 
  for(i=nfreeze;i<n;++i) rn += m[i];
  rn = 1./rn;
  for(k=0;k<DIM;++k) tv[k]=0.;
  for(i=nfreeze;i<n;++i) for(k=0;k<DIM;k++) tv[k] += v[i*DIM+k]*m[i];
  for(k=0;k<DIM;++k) tv[k] *= rn;
  for(i=nfreeze;i<n;++i) for(k=0;k<DIM;k++) v[i*DIM+k] -= tv[k];
}
