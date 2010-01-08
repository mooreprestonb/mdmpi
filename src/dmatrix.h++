
#ifndef _DMATRIX_H
#define _DMATRIX_H
int *ivector(int nl,int nh);
void free_dvector(double *v,int nl,int nh);
double *dvector(int nl,int nh);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
float **fmatrix(int nrl,int nrh,int ncl,int nch);
void free_fmatrix(float **m,int nrl,int nrh,int ncl,int nch);
double ***d3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh);
#endif
