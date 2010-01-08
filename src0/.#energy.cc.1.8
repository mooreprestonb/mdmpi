/*! \file energy.cc
  file which holds energy routines
*/

//! prototype to write to errorlog 
extern void writeErr(const char *s); 

#include "energy.h++"

ENERGY::ENERGY(void) //! constructor
{
  _dof=_navg=0;
  _peng=_keng=_spe=_ske=0.0;
  _plr = _psr = 0.;
  _lj=_elsr=_ellr=_elcorr=0;
}

double ENERGY::addAvg(double ke,double pe,double sr,double lr,
		      double plj,double esr,double elr,double ecorr)
{
  _navg++;
  _peng = pe;  _keng = ke; 
  _spe += pe;  _ske += ke;
  _psr = sr; _plr = lr;
  _lj = plj; _elsr = esr; _ellr = elr; _elcorr = ecorr;
  return pe+ke;
}

void ENERGY::printavg(FILE *fp)
{
  // #define OLD
#ifdef OLD
  fprintf(fp,"<HAM>= % 10.5e ~ % 10.5e\n",avgHAM(),ham());
  fprintf(fp,"<KE> = % 10.5e ~ % 10.5e\n",avgKE(),ke());
  fprintf(fp,"<PE> = % 10.5e ~ % 10.5e\n",avgPE(),pe());
#else
  fprintf(fp,"<energyAVG>\n");
  fprintf(fp,"<ham> % 10.5e ~ % 10.5e </ham>\n",avgHAM(),ham());
  fprintf(fp,"<ke> % 10.5e ~ % 10.5e </ke> \n",avgKE(),ke());
  fprintf(fp,"<pe> % 10.5e ~ % 10.5e </pe>\n",avgPE(),pe());
  fprintf(fp,"</energyAVG>\n");
#endif
}

void ENERGY::print(FILE *fp)
{
  fprintf(fp,"<energy> <dof> %d </dof> <navg> %d </navg> ",_dof,_navg);
  fprintf(fp,"<tke> %g </tke> <tpe> %g </tpe> </energy>",_ske,_spe);
  // fprintf(fp,"\tke=%g, pe=%g, lj=%g, elsr=%g, ellr=%g, elcorr=%g\n",
  // keng,peng,_lj,_elsr,_ellr,_elcorr);
}

void ENERGY::iprint(int istep,FILE *fp)
{
  double rd = 1./_dof;
  double ham = _keng+_peng;
  double elec = _elsr+_ellr+_elcorr;
  // rd = 1;
  fprintf(fp,"%d %g %g %g %g %g %g %g %g %g %g %g\n",istep,ham*rd,
	  _keng*rd,_peng*rd,_psr*rd,_plr*rd,
	  _lj*rd,elec*rd,_elsr*rd,_ellr*rd,_elcorr*rd,rd);
}

void ENERGY::read(int nmol,const std::string & name)
{
  int i;
  const int MAXLINE=256;
  char line[MAXLINE];
  std::string sline;
  FILE *fp;
  if((fp = fopen(name.c_str(),"r"))==NULL){
    sprintf(line,"ERROR: can't open coordfile %s",name.c_str());
    sline = line;writeErr(sline.c_str()); exit(1);
  }
  if(fgets(line,MAXLINE,fp)==NULL) { // header
    sprintf(line,"ERROR: empty file %s?",name.c_str());
    sline = line;writeErr(line); exit(1);
  }
  for(i=0;i<nmol;i++){
    if(fgets(line,MAXLINE,fp)==NULL) { // atoms
      sprintf(line,"ERROR: in coordfile file %s? at line %d",name.c_str(),i+1);
      sline = line;writeErr(line); exit(1);
    }
  }
  if(fgets(line,MAXLINE,fp)==NULL) { // cell 
    sprintf(line,"ERROR: can't read cell in fine %s?",name.c_str());
    sline = line;writeErr(line); exit(1);
  }
  if(fgets(line,MAXLINE,fp)==NULL) { // energy
    sprintf(line,"ERROR: in file %s? while reading in Energy\n",name.c_str());
    sline = line;writeErr(line); 
    sprintf(line,"Mayby an initial file instead of a restart file?");
    sline = line;writeErr(line);     exit(1);
  } else {
    sscanf(line,"Energy: dof= %d, navg= %d, tke= %lg tpe = %lg\n",
	  &_dof,&_navg,&_ske,&_spe);
  }
  fclose(fp);
}
