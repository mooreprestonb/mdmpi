/*! \file mdmol.cc
  Main C++ program to do Molecular Dynamics 
*/

#include "mdmol.h++"
#include "flags.H"

/*! \fn freeScratch()
  free the scratch arrays that have been used
*/
void freeScratch(void){
  freeLJScr();
  freeEwaldStatic();
  // freeEnergyScr();
}

void propeller(void)
{
  static int iwrite=0;
  fprintf(stdout,"\r%75s\r"," ");
  switch(iwrite++%4){
  case 0: fprintf(stdout,"\\"); break;case 1: fprintf(stdout,"|"); break;
  case 2: fprintf(stdout,"/"); break; case 3: fprintf(stdout,"-"); break;
  }
  fprintf(stdout," ");
}

void progress(int istep, int nstep,double wtime)
{
  double ris,tfinish; 
  propeller();
  if(istep==0) {  ris = 1.;   tfinish = 0;} 
  else { ris = 1./double(istep);  tfinish = wtime*nstep/double(istep); }
  fprintf(stdout,"%d %4.0f%% cpu/step= %6.4g, tot=%6g sec (~%6g tot %6g left)",
	  istep,100.*istep/double(nstep),wtime*ris,wtime,tfinish,tfinish-wtime); 
  fflush(stdout);
}

/*! \fn main(int argc,char *argv[])
  main routine that does everything
*/
int main(int argc,char *argv[])
{
  int istep;
  SIMVARS simvars;
  ENERGY energy;
  COORDS coords;
  NEIGHBOR ngbr;
  TIMING time;
  FLAGS flags;

  flags.read(argc,argv);
  setup(flags,simvars,coords,ngbr,energy);
  time.start();

  istep = simvars.istep();
  while(istep<=simvars.nstep()){
    coords.reset();
    average(simvars,coords,ngbr,energy);
    if(simvars.rank()==0){
      if((istep%simvars.iwrite())==0){
	if(!simvars.hamname().empty())
	  storeHam(simvars.hamname().c_str(),istep,energy,simvars.rank());
	if(!simvars.restname().empty()) saveRestart(simvars,coords,ngbr,energy);
	if(!simvars.trajname().empty()) storeConf(simvars,coords,ngbr);
	writeELog(istep,simvars,energy);
	if(!flags.q()) progress(istep,simvars.nstep(),time.totalwalltime());
      }
    }
    stepvv(simvars,coords,ngbr);      /* integrate VV*/
    simvars.istep(++istep);
  }
  saveRestart(simvars,coords,ngbr,energy);
  char line[MAXLINE];
  sprintf(line,"\nTiming of %s on %s (id=%d)",argv[0],simvars.hostname().c_str(),simvars.pid()); 
  writeOut(line);
  sprintf(line,"(id=%d) cputime = %g sec, walltime = %g sec",simvars.pid(),time.cputime(),time.walltime()); 
  writeOut(line);

  freeScratch();
  if(!flags.q()) fprintf(stdout,"\n");
  return 0;
}

