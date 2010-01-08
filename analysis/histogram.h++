
#ifndef _HISTOGRAM2D_H
#define _HISTOGRAM_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

class HISTOGRAM
{
private:
  int _npt,_nadd,_ibin,_nmax;
  char *_filename;
  double _dx,_rmax,_rmin,_xmin,_xmax,*_hist;
public:
  HISTOGRAM(void);
  ~HISTOGRAM(void);
  char * filename(void){return _filename;}
  double dx(void){return _dx;}
  double* hist(void){return _hist;}
  int ntp(void){return _npt;}
  void set(int,double,double,char *);
  void resize(int);
  void bin(double);
  void print(void);
  double xmax(void){return _xmax;}
  double xmin(void){return _xmin;}
};

#endif
