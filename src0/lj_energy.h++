
#ifndef _LJ_ENERGY_H
#define _LJ_ENERGY_H

double ljspot(void);
void freeLJScr(void);
void LJsetup(int,double[],double[],double,double);
double pot_lj(int,double[],int[]);
double force_lj(int,double[],int[]);

#endif
