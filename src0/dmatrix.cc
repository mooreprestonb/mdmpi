/* subroutines to return memory of matrix or tensor3 */

#include <stdlib.h>
#include <stdio.h>

#include "dmatrix.h++"

/*---------------------------------------------------------------------*/
/* allocate an integer vector with subscipt range v[nl...nh] */
int *ivector(int nl,int nh)
{
  int *v; 
  v = (int *)malloc((nh-nl+1)*sizeof(int));
  return v-nl;
}
void free_ivector(int *v,int nl,int nh){ free(v+nl);}
/*---------------------------------------------------------------------*/
/* allocate an double vector with subscipt range v[nl...nh] */
double *dvector(int nl,int nh)
{
  double *v; 
  v = (double *)malloc((nh-nl+1)*sizeof(double));
  return v-nl;
}
void free_dvector(double *v,int nl,int nh){ free(v+nl);}
/*---------------------------------------------------------------------*/
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(int nrl,int nrh,int ncl,int nch) {
  int i,nrow,ncol;
  double **m;

  nrow = nrh-nrl+1;  ncol = nch-ncl+1;

  if(nrow*ncol==0) return (double **)NULL;
  /* allocate pointers to rows */

  m = (double **)malloc(nrow*sizeof(double*));
  if(!m){
    fprintf(stderr,"ERROR: can't allocate memory for matrix\n");
    exit(1);}
  m -= nrl;

  /* allocate rows */
  m[nrl] = (double *)calloc(nrow*ncol,sizeof(double));
  if(!m[nrl]){
    fprintf(stderr,"ERROR: can't allocate memory for arrays for matrix\n");
    exit(1);
  }
  m[nrl] -= ncl;

  /* set pointers to rows */
  for(i=nrl+1;i<=nrh;i++) m[i] = m[i-1]+ncol;

  return m;
}
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch){
  free((char *)(m[nrl]+ncl)); free((char *)(m+nrl));
}

/*----------------------------------------------------------------------*/
float **fmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow,ncol;
  float **m;

  nrow = nrh-nrl+1;  ncol = nch-ncl+1;

  if(nrow*ncol==0) return (float **)NULL;
  /* allocate pointers to rows */

  m = (float **)malloc(nrow*sizeof(float*));
  if(!m){
    fprintf(stderr,"ERROR: can't allocate memory for matrix\n");
    exit(1);}
  m -= nrl;

  /* allocate rows */
  m[nrl] = (float *)calloc(nrow*ncol,sizeof(float));
  if(!m[nrl]){
    fprintf(stderr,"ERROR: can't allocate memory for arrays for matrix\n");
    exit(1);
  }
  m[nrl] -= ncl;

  /* set pointers to rows */
  for(i=nrl+1;i<=nrh;i++) m[i] = m[i-1]+ncol;

  return m;
}

/*----------------------------------------------------------------------*/
void free_fmatrix(float **m,int nrl,int nrh,int ncl,int nch)
{
  free((char *)(m[nrl]+ncl)); free((char *)(m+nrl));
}

/*----------------------------------------------------------------------*/

double ***d3tensor(int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
  int i,j,nrow,ncol,ndep;
  double ***t;

  nrow = nrh-nrl+1; ncol = nch-ncl+1; ndep = ndh-ndl+1;

  /* allocate top pointer */

  t = (double ***)malloc(nrow*sizeof(double**));
  if(!t){
    fprintf(stderr,"ERROR: (1) can't allocate memory for tensor\n");
    exit(1);
  }
  t -= nrl;

  /* allocate matrix pointers */
  t[nrl] = (double **)malloc(nrow*ncol*sizeof(double*));
  if(!t[nrl]){
    fprintf(stderr,"ERROR: (2) can't allocate memory for arrays in tensor\n");
    exit(1);
  }
  t[nrl] -= ncl;

  /* allocate final pointers and set all pointers */
  t[nrl][ncl] = (double *)calloc(nrow*ncol*ndep,sizeof(double));
  if(!t[nrl][ncl]){
    fprintf(stderr,"ERROR: (3) can't allocate memory for rows in tensor\n");
    exit(1);
  }
  t[nrl][ncl] -= ndl;


  /* set pointers to rows */
  for(j=ncl+1;j<=nch;j++) t[nrl][j] = t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++){
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j] = t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */

  return t;
}
/*----------------------------------------------------------------------*/

void free_d3tensor(double ***t,int nrl,int nrh,int ncl,int nch,int ndl,int ndh)
{
  free((char *)(t[nrl][ncl]+ndl)); 
  free((char *)(t[nrl]+ncl));
  free((char *)(t+nrl));
}
/*----------------------------------------------------------------------*/
