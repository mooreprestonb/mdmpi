
#include <math.h>
#include "mathutil.h++"

double dot2(int n,double a[])
{
  int i;
  double sum=0.;
  for(i=0;i<n;i++) sum += a[i]*a[i];
  return sum;
}

double dot(int n,double a[],double b[])
{
  int i;
  double sum=0.;
  for(i=0;i<n;i++) sum += a[i]*b[i];
  return sum;
}
