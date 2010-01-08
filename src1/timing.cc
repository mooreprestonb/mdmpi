/* class for time processes */

#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

#define MAXTICS 2147483647 // int max ;-)

#ifdef PARA
#include <mpi.h>
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC  1000000
#endif

#include "timing.h++"

/*-------------------------------------------------*/
void TIMING :: start(void)
{
  cticsBegin = cticsOld = cticsNew = clock();
  ttime = 0; itime=0;
#ifdef PARA
  secBegin = secOld = secNew = double(MPI_Wtime());
#else 
  secBegin = secOld = secNew = double(time(&ttime));
#endif
}

/*-------------------------------------------------*/
TIMING :: TIMING(void){start();}
/*-------------------------------------------------*/
TIMING :: ~TIMING(void){}
/*-------------------------------------------------*/
double TIMING :: cputime(void)
{
  
  cticsNew = clock();
  long ctics; // ctics could rap around so we max out at ~35min (2147sec).
  if(cticsNew<cticsOld) ctics = (MAXTICS - cticsOld)+cticsNew;
  else ctics = cticsNew-cticsOld;
  cticsOld = cticsNew;
  return (double(ctics)/double(CLOCKS_PER_SEC));
}

double TIMING :: totalcputime(void)
{
  long cticsNow = clock();
  long ctics; // ctics could rap around so we max out at ~35min (2147sec).
  if(cticsNow<cticsBegin) ctics = (MAXTICS - cticsBegin)+cticsNew;
  else ctics = cticsNow-cticsBegin;
  return (double(ctics)/double(CLOCKS_PER_SEC));
}
/*-------------------------------------------------*/
double TIMING :: walltime()
{
#ifdef PARA
  secNew = double(MPI_Wtime());
#else 
  secNew = double(time(&ttime));
#endif
  double dtime = secNew-secOld;
  secOld = secNew;
  return dtime;
}
/*-------------------------------------------------*/
double TIMING :: totalwalltime()
{
#ifdef PARA
  secNew = double(MPI_Wtime());
#else 
  secNew = double(time(&ttime));
#endif
  double dtime = secNew-secBegin;
  return dtime;
}
/*-------------------------------------------------*/
