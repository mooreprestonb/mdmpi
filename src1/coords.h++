/*! \file coords.h++
  input file and definition of the coords class
*/

#ifndef _COORDS_H
#define _COORDS_H

#include <stdlib.h>
#include <stdio.h>
#include <string>

void writeOut(const char *s);
void writeLog(const char *s);

#include "mdinclude.h++"

//! the coords class
/*!
  This class hold the coordinates for the simulation <br> 
  (positions, velocities, accelerations, charges, and masses)
*/

class COORDS
{
private:
  int _natoms,_ntypes;
  int *_itype;
  int *_icells; //!< number of cells each coord had moved
  double _box[DIM];  // lenght of box;
  double *_x,*_v,*_a,*_qch,*_mass;
public:
  COORDS(void);   //!< constructor
  ~COORDS(void);  //!< destructor
  //! return the number of atoms
  int natoms(void){return _natoms;}
  //! set and return the number of atoms
  int natoms(int n){_natoms=n;return _natoms;}
  //! return the number of types
  int ntypes(void){return _ntypes;}
  //! set and return the number of types
  int ntypes(int n){_ntypes=n;return _ntypes;}
  //! return pointer to positions 
  double* x(void) {return _x;};
  //! return pointer to velocities
  double* v(void) {return _v;};
  //! return pointer to accelerations 
  double* a(void) {return _a;};
  //! return pointer to charges
  double* qch(void) {return _qch;};
  //! return pointer to masses
  double* mass(void) {return _mass;};
  //! return index of atom types 
  int * itype(void) {return _itype;};
  //! read pos & vel from file
  int read(int,const std::string &);
  //! prints positions & velocities
  void print(FILE * = stdout);
  //! pull atoms back into the box (do we want to scale as well?)
  void reset(void);
};

#endif
