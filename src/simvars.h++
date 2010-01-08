/*! \file simvars.h++
   file that defines the simvars class <br>
   this class holds all the simulation varibles (ie simvars)
   for running the molecular dynamics 
   \class simvars
*/

#ifndef _SIMVARS_H
#define _SIMVARS_H

#include <string>

class SIMVARS 
{
private:
  std::string _hostname,_inname,_restname;
  std::string _hamname,_coordname,_trajname;
  std::string _logname,_errname;
  std::string _paramfile,_systemfile;
  int _nmol,_istep,_nstep,_ipair,_iperd;
  int _ncell,_rank,_size,_dof,_pid,_kmax,_nfreeze;
  double _dt,_xbox,_rcut,_skin,_alpha,_esurf;
  int _zerotm,_resetpos,_writebox,_cutoff;
  int _iwrite,_calctype;
public:
  SIMVARS(void); // constructor
  ~SIMVARS(void); // destructor
  //! set hostname
  void hostname(const std::string &t){_hostname = t;}  
  //! return hostname
  std::string hostname(void){return _hostname;}
  //! set input file name
  void inname(const std::string &t){_inname = t;}
  //! return input file name
  std::string inname(void){return _inname;}
  //! set restart file name
  void restname(const std::string &t){_restname = t;}
  //! return restart file name
  std::string restname(void){return _restname;}
  //! set trajectory file name 
  void trajname(const std::string &t){_trajname = t;}
  //! return trajectory name 
  std::string trajname(void){return _trajname;}
  //! set coordinate file name (file that hold initial coords)
  void coordname(const std::string &t){_coordname = t;} 
  //! return coordinate file name (file that hold initial coords)
  std::string coordname(void){return _coordname;}
  //! set the hamiltonian name (file that hold the energies)
  void hamname(const std::string &t){_hamname = t;} 
  //! return the hamiltonian name (file that hold the energies)
  std::string hamname(void){return _hamname;}
  //! set the name of the log file
  void logname(const std::string &t){_logname = t;}   
  //! return the name of the log file
  std::string logname(void){return _logname;}
  //! set the name of the error file
  void errname(const std::string &t){_errname = t;}   
  //! return the name of the error file
  std::string errname(void){return _errname;}
  //! set the name of the parameter file
  void paramfile(const std::string &t){_paramfile = t;}   
  //! return the name of the param file
  std::string paramfile(void){return _paramfile;}
  //! set the name of the error file
  void systemfile(const std::string &t){_systemfile = t;}   
  //! return the name of the system file
  std::string systemfile(void){return _systemfile;}
  //! set the process id
  void pid(int p){_pid=p;} 
  //! return the process id
  int pid(void){return _pid;}
  //! return the number of atoms
  int nmol(void){return _nmol;}
  //! return the number of total steps to be taken 
  int nstep(void){return _nstep;}
  //! return the number of steps we are on
  int istep(void){return _istep;} 
  //! return/set the number of steps we are on
  int istep(int is){_istep=is;return _istep;}
  //! return the neighbor pair type
  int ipair(void){return _ipair;}
  //! return the periodicity
  int iperd(void){return _iperd;}
  //! return the number of link list cells we are using
  int ncell(void){return _ncell;}
  //! return the number of frozen atoms
  int nfreeze(void){return _nfreeze;}
  //! return the processor rank
  int rank(void){return _rank;}
  //! return the processor rank
  int rank(int rank){_rank = rank;return _rank;}
  //! return the processors size
  int size(void){return _size;}
  //! return the degrees of freedom
  int size(int size){_size=size;return _size;}
  //! return the degrees of freedom
  int dof(void){return _dof;}
  //! return kmax
  int kmax(void){return _kmax;}
  //! return whether or not to zero the total momentum
  int zerotm(void){return _zerotm;} 
  //! return the cutoff value
  int cutoff(void){return _cutoff;} 
  //! return if we want to reset the positions on writing out
  int resetpos(void){return _resetpos;}
  //! return if we want to write the box or not
  int writebox(void) {return _writebox;}
  //! return what type of calculation we are doing (initial or continue)
  int calctype(void){return _calctype;}
  //! off often to write out
  int iwrite(void){return _iwrite;}
  // int iwrite(int iw){_iwrite=iw;return _iwrite;}
  //! the cut of distance
  double rcut(void){return _rcut;}
  //! the skin distance
  double skin(void){return _skin;}
  //! the alpha
  double alpha(void){return _alpha;}
  //! the electro static purativitiy of the surrounding medic, us in ewald
  double esurf(void){return _esurf;}
  //! the time step
  double dt(void){return _dt;}
  //! the box length
  double xbox(void){return _xbox;}
  //! read the input simulation parameters
  void readinit(void);
  //! print out the class variables
  void print(FILE * = stdout);
};

#endif
