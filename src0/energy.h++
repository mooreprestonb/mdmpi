/*! \file energy.h++
  file with functions for the energy class
*/

#ifndef _ENERGY_H
#define _ENERGY_H

#include <stdlib.h>
#include <stdio.h>
#include <string>

//! the energy class
/*!
  This class hold the varibles and sums the different energy contributions 
  <br> like the kinetic and potential, electrostatic etc...
*/

class ENERGY
{
private:
  int _dof,_navg;
  double _peng,_keng,_ske,_spe,_plr,_psr;
  double _lj,_elsr,_ellr,_elcorr;
public:
  ENERGY(void);
  // ~ENERGY(void);
  //! return the number of degrees of freedom ie.. natoms*Dim - nfreeze*DIM
  int dof(void){return _dof;}
  //! set the number of degrees
  int dof(int df){_dof = df; return _dof;}
  //! return the keng
  double ke(void){return _keng;}
  //! store and return the keng
  double ke(double eng){_keng = eng; return _keng;}
  //! return the potential energy
  double pe(void){return _peng;}
  //! store potential energy
  double pe(double tpe){_peng = tpe;  return _peng;}
  //! return the total energy
  double ham(void){return (_keng+_peng);}
  //! return the average kinetic energy
  double avgKE(void){return _ske/_navg;}
  //! return the average potential energy
  double avgPE(void){return _spe/_navg;}
  //! return the average total energy
  double avgHAM(void){return (_ske+_spe)/_navg;}
  //! add to the energy variable
  double addAvg(double,double,double,double,double,double,double,double);
  //! print the averge energy to fp;
  void printavg(FILE *fp=stdout);
  //! print the energy class out;
  void print(FILE *fp=stdout);
  //! print the instantanous values of the energy class to fp;
  void iprint(int, FILE *fp=stdout);
  //! read int the energy for restart file;
  void read(int nmol,const std::string &);
};

#endif
