/*! \file simvars.cc 
  functions for the simvars class
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #define DEBUG
const int DIM=3; //!< Dimensionality of the system

#include <list>
#include <algorithm>

#include "simkeys.h++"
#include "simvars.h++"

SIMVARS :: ~SIMVARS(void){}

SIMVARS :: SIMVARS(void)
{
  _rank=0; _size=1; _pid=-1;
  _nmol=0;_nstep=0;_istep=0;
  _ipair=2;_iperd=0;_ncell=3;_dof=-1;
  _kmax=10;_nfreeze=0;
  _dt=1.;_xbox=1.;
  _rcut=3.0;_skin=1.0,_alpha=2.;
  _zerotm = _resetpos = _writebox = 1;// on by default
  _cutoff = 0; // off by default
  _iwrite = _calctype = 0;
  _esurf = -1; // use -1 as infinity of ewald surface dielectic constant
  // now strings so we don't need to initialize
  //_hostname = _inname = _restname = _hamname = (char *)NULL;
  //_trajname = _coordname = _logname = _errname = (char *)NULL;
}

/*! \fn void SIMVARS :: print(FILE *fp)
  print all the external variable of the SIMVARS class
  \param fp file pointer to print to
*/

static void printname(FILE *fp,const char *name,const char * _name)
{
  if(_name!=NULL) fprintf(fp,"%s = \"%s\"",name,_name);
  else fprintf(fp,"%s = \"\"",name);
}


void SIMVARS :: print(FILE *fp)
{
  fprintf(fp,"nmol = %d, istep = %d, nstep = %d, iperd=%d, xbox=%g, dt=%g\n",
	  _nmol,_istep,_nstep,_iperd,_xbox,_dt);
  fprintf(fp,"rank = %d, size = %d, write = %d, itype = %d\n",
	  _rank,_size,_iwrite,_calctype);
  fprintf(fp,"ipair=%d, ncell=%d, rcut = %g, skin = %g\n",
	  _ipair,_ncell,_rcut,_skin);
  fprintf(fp,"kmax = %d, alpha = %g, esurf = %g, icutoff = %d, nfreeze = %d\n",
	  _kmax,_alpha,_esurf,_cutoff,_nfreeze);
  fprintf(fp,"hostname = %s, pid = %d\n",_hostname.c_str(),_pid);
  fprintf(fp,"zerotm = %d, resetpos = %d, writebox = %d\n",
	  _zerotm,_resetpos,_writebox);
  printname(fp,"setname",_inname.c_str()); fprintf(fp,", ");
  printname(fp,"restname",_restname.c_str()); fprintf(fp,"\n");
  printname(fp,"coordname",_coordname.c_str()); fprintf(fp,", ");
  printname(fp,"trajname",_trajname.c_str()); fprintf(fp,", ");
  printname(fp,"hamname",_hamname.c_str());  fprintf(fp,"\n");
  printname(fp,"logname",_logname.c_str());fprintf(fp,", ");
  printname(fp,"errname",_errname.c_str()); fprintf(fp,"\n");
  printname(fp,"systemfile",_systemfile.c_str()); fprintf(fp,"\n");
  printname(fp,"paramfile",_paramfile.c_str()); fprintf(fp,"\n");
}


