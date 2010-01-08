/*! /file energy_intra.cc
  file with routines to calculate intramolecular 
  energies/forces/hessian we need
  bonds bends torsions and cross terms if needed
*/

#include "mdmol.h++"

static int nxdis=0;
static double *xdis=NULL;
