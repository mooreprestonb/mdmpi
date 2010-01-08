
#ifndef _TIMING_H
#define _TIMING_H

#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

class TIMING
{
private:
  int cticsOld,cticsNew,cticsBegin;
  double secOld,secNew,secBegin;
  clock_t itime;
  time_t  ttime;
public:
  TIMING(void);
  ~TIMING(void);
  void start(void);
  double cputime(void);
  double walltime(void);
  double totalwalltime(void);
  double totalcputime(void);
};

#endif
