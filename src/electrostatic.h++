/*! \file electrostatic.h++
  file to hold the prototypes and meta includes for 
  the subroutines for electrostatic interactions
  both ewald and real space terms.
*/

#ifndef _ELECTROSTATIC_H
#define _ELECTROSTATIC_H

double eleclrpot(void);
double elecsrpot(void);
double eleccorrpot(void);
double potElec(int,double[],double[]);
double forceElec(int,double[],double[]);
double k_ewald(int,double [],double []);
double fk_ewald(int,double [],double [],double []);
void freeEwaldStatic(void);
void ewaldSetup(int,double *,int,double,double,double,int,int);
void printEcorr(void);
void addEcorr(int,int);
double ecorr(double *px,double *fx,double *qch);
double potEcorr(double *px,double *qch);
void cellDipole(int n,double *x,double *q,double d[3]);
double cellQuadrupole(int n,double *x,double *q);
void potElecx(int nch,double *xi,double *qi, double x[]);
double directSum(const int nch,const double px[],const double qch[]);
#endif
