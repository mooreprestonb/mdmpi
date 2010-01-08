#ifndef MAX
#define MAX(A,B) ((A>B)?(A):(B))
#endif
#ifndef MIN
#define MIN(A,B) ((A<B)?(A):(B))
#endif

#include "histogram.h++"

HISTOGRAM :: HISTOGRAM(void)
{
  _npt = 2;  _nadd = 0; _ibin=1; _nmax=0;
  _rmin = _rmax = 0.;
  _xmin = _xmax = 0.;
  _dx = 0;
  _hist = new double[_npt];
  for(int i=0;i<_npt;i++) _hist[i] = 0;
  _filename = NULL;
}

HISTOGRAM :: ~HISTOGRAM(void)
{
  delete[] _hist;
  if(_filename != NULL) free(_filename);
}

void HISTOGRAM::resize(int ns)
{
  if(ns<2) ns=2;
  if(_npt != ns) delete[] _hist;
  _npt = ns;
  _hist = new double[_npt];
  for(int i=0;i<_npt;i++) _hist[i] = 0;
  if(_ibin) _nmax=_npt-2;
  else _nmax=_npt-1;
}

void HISTOGRAM :: set(int size,double xmin,double xmax,char *name)
{
  if(size<2) size=2;
  resize(size);
  _nadd = 0;
  _xmin = xmin;  _xmax = xmax;
  _dx = (_xmax-_xmin)/(double)(_npt-1);
  for(int i=0;i<_npt;i++) _hist[i] = 0;
  if(name!=NULL) _filename = strdup(name);
}

void HISTOGRAM::bin(double pt)
{
  int n;
  if(_nadd++ == 0) _rmax=_rmin=pt;

  _rmax = MAX(_rmax,pt);
  _rmin = MIN(_rmin,pt);
  pt = (pt-_xmin)/_dx;
  n = (int)pt;    n = MAX(n,0);    n = MIN(n,_nmax);
  if(_ibin){ 
    pt -= n;       pt = MAX(pt,0.); pt = MIN(pt,1.);
    _hist[n] += 1. - pt; _hist[n+1] += pt;
  } else _hist[n]++;
}

void HISTOGRAM::print(void)
{
  int i;
  FILE *fp = fopen(_filename,"w");

  if(fp == NULL){
    fprintf(stderr,"ERROR: can't open %s\n",_filename);
    return;
  }
  fprintf(fp,"# npt=%d, min=%g, max=%g, n=%d, rmin=%g, rmax=%g\n",
	  _npt,_xmin,_xmax,_nadd,_rmin,_rmax);
  double pt = 1./(_dx*(double)(_nadd));
  for(i=0;i<_npt;i++){
    fprintf(fp,"%g %g\n",i*_dx+_xmin,_hist[i]*pt);
  }
  fclose(fp);
}

