/*! \file mdOut.cc
  file to hold the input output of mdmol
*/

#include "mdmol.h++"

static FILE *fpl=stdout;  //!< logfile
static FILE *fpo=stdout;  //!< output
static FILE *fpe=stderr;  //!< errfile

//! set stdout pointer
void setStdout(const char *file){
  if(fpo != stdout ){fprintf(fpo,"</outfile>\n"); fclose(fpo);}
  if(file == NULL || strcmp(file,"")==0) fpo = stdout;
  else {
    if((fpo = fopen(file,"w"))==NULL){
      fprintf(fpe,"Can't set stdout to \"%s\", setting to stdout\n",file);
      fpo = stdout;
    }
    fprintf(fpo,"<outfile>\n");
  }
}

//! set stderr pointer
void setStderr(const char *file){
  if(fpe != stderr ){ fprintf(fpe,"</errorfile>\n"); fclose(fpe);}
  if(file == NULL || strcmp(file,"")==0) fpe = stderr;
  else {
    if((fpe = fopen(file,"w"))==NULL){
      fprintf(stderr,"Can't set stderr to \"%s\", setting to stderr\n",file);
      fpe = stderr;
    }
    fprintf(fpe,"<errorfile>\n");
  }
}

//! set logfile file pointer
FILE * setLogfile(const char *file){
  if(fpl != stdout ){ 
    fprintf(fpo,"</logfile>\n");
    fclose(fpo);
  }
  if(file == NULL || strcmp(file,"")==0) fpl = stdout;
  else
    if((fpl = fopen(file,"w"))==NULL){
      fprintf(stderr,"Can't set logfile to \"%s\", setting to stdout\n",file);
      fpl = stdout;
    }
  fpo = fpl; // set output to log
  // fprintf(fpo,"<logfile>\n");
  return fpl;
}

//! write to stdout
void writeOut(const char *s){fprintf(fpo,"%s\n",s); fflush(fpo);}
//! write to stderr
void writeErr(const char *s){fprintf(fpe,"%s\n",s); fflush(fpe);}
//! write to log
void writeLog(const char *s){fprintf(fpl,"%s\n",s); fflush(fpl);}

//! write energy log
void writeELog(int istep,SIMVARS &simvars,ENERGY &energy)
{
  // printf("%d %g %g %g\n",istep,energy.ham(),energy.ke(),energy.pe());
  if(fpl!=NULL){
    // fprintf(fp,"\r%75s\r%8d "," ",istep);
    fprintf(fpl,"<step> %8d \n",istep);
    //fprintf(fpl,"\nStep = %8d\n",istep);
    energy.printavg(fpl);
    fprintf(fpl,"</step>\n");
    fflush(fpl);
  }
}

void saveRestart(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr,ENERGY &energy)
{
  if(simvars.restname().empty()){
    fprintf(fpe,"Warning: restart file name is null, not saving!\n");
    return;
  }

  double *x = coords.x();
  double *v = coords.v();    
  double *qch = coords.qch();
  double *m = coords.mass();
  ngbr.getcm(x,v);

  if(simvars.rank()==0){
    int i,j,k;
    FILE *fpr;
    int nmol = simvars.nmol();
    double *xt = (double *)malloc(2*nmol*DIM*sizeof(double));
    double *vt = xt + nmol*DIM;
    for(i=0;i<nmol*DIM;i++) {xt[i] = x[i];vt[i] = v[i];}
    if((fpr=fopen(simvars.restname().c_str(),"w"))==NULL){
      fprintf(fpe,"ERROR: can't open %s for writing restart\n",
	      simvars.restname().c_str());
      exit(1);
    }
    if(simvars.zerotm()) zerotm(nmol,DIM,simvars.nfreeze(),vt,m);
    if(simvars.resetpos()) periodic(nmol,xt);

    fprintf(fpr,"<coords>\n");
    fprintf(fpr,"<natoms> %d </natoms>\n",nmol);
    fprintf(fpr,"<step> %d </step>\n", simvars.istep());
    fprintf(fpr,"<timestep> %g </timestep>\n", simvars.dt());
    fprintf(fpr,"<boxsize> %g </boxsize>\n",simvars.xbox());
    fprintf(fpr,"<coordFileName> %s </coordFileName>\n",
	    simvars.coordname().c_str());
    fprintf(fpr,"<restartFileName> %s</restartFileName>\n",
	    simvars.restname().c_str());
    fprintf(fpr,"<energyFileName> %s </energyFileName>\n",
	    simvars.hamname().c_str());
    fprintf(fpr,"<trajectorFileName> %s </trajectorFileName>\n",
	    simvars.trajname().c_str());

    for(i=0;i<nmol;i++){
      j = i*DIM;
      fprintf(fpr,"<xyz id=\"%d\"> ",i);
      for(k=0;k<DIM;k++) fprintf(fpr,"%g ",xt[j+k]); fprintf(fpr,"\t");
      fprintf(fpr,"</xyz> \n <vxyz id=\"%d\"> ",i);
      for(k=0;k<DIM;k++) fprintf(fpr,"%g ",vt[j+k]); fprintf(fpr,"\t");
      fprintf(fpr,"</vxyz> \n <charge id=\"%d\"> ",i);
      fprintf(fpr,"%g </charge> <mass id=\"%d\"> %g </mass>\n",qch[i],i,m[i]);
    }
    free(xt);
    if(simvars.writebox()) 
      fprintf(fpr,"<hmatrix> %g 0 0  0 %g 0  0 0 %g </hmatrix>\n",
	      simvars.xbox(),simvars.xbox(),simvars.xbox());
    energy.print(fpr);
    fprintf(fpr,"</coords>\n");
    fclose(fpr);
  }
}

void storeConf(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr)
{
  static int istart=0;
  static FILE *fpt=NULL;
  int i,j,k,nmol;

  double *x = coords.x();
  double *v = coords.v();
  double *xt;
  ngbr.getcm(x,v);

  if(simvars.rank()==0){
    nmol = simvars.nmol();
  
    if(simvars.trajname().empty()){
      fprintf(fpe,"Warning: trajectory file name is null, not saving!\n");
      return;
    }
    if(istart==0){
      if((fpt = fopen(simvars.trajname().c_str(),"w"))==NULL){
	fprintf(fpe,"ERROR: can't open trajectory file \"%s\" \n",
		 simvars.trajname().c_str());
	exit(1);
      }
      istart=1;
      fprintf(fpt,"<trajectoryFile>\n");
    }

    fprintf(fpt,"<coords>\n");
    fprintf(fpt,"<nmol> %d </nmol> <step> %d </step> <time> %g </time>\n",
	    nmol,simvars.istep(),simvars.istep()*simvars.dt());

    if(simvars.resetpos()) {
      xt = (double *)malloc(nmol*DIM*sizeof(double));
      for(i=0;i<nmol*DIM;i++) {xt[i] = x[i];}
      periodic(nmol,xt);
    } else {
      xt = x;
    }

    for(j=i=0;i<nmol;i++,j+=DIM) {
      fprintf(fpt,"<xyz id=\"%d\"> ",i);
      for(k=0;k<DIM;k++) fprintf(fpt,"%g ",xt[j+k]);
      fprintf(fpt,"</xyz> ");
      fprintf(fpt,"<vxyz id=\"%d\"> ",i);
      for(k=0;k<DIM;k++) fprintf(fpt,"%g ",v[j+k]);
      fprintf(fpt,"</vxyz>\n");
    }

    if(simvars.resetpos()) {
      free(xt);
    }

    if(simvars.writebox()){
      fprintf(fpt,"<hmatrix> %g 0 0  0 %g 0  0 0 %g </hmatrix>\n",
	      simvars.xbox(),simvars.xbox(),simvars.xbox());
    }
    fprintf(fpt,"</coords>\n");
  }
}

void storeHam(const char *filen,int istep,ENERGY &energy,int rank)
{
  static int istart=0;
  static FILE *fph;
  
  if(rank==0){
    if(filen==NULL){
      fprintf(fpe,"Warning: energy file name is NULL, not saving!\n");
      return;
    }
    if(istart==0){
      istart=1;
      if((fph = fopen(filen,"w"))==NULL){
	fprintf(fpe,"ERROR: can't open energy file\"%s\", not saving\n",filen);
	// exit(1);
      }
      if(fph != NULL) {
	fprintf(fph,"# istep d=%d, ham ke pe psr plr ",energy.dof());
	fprintf(fph,"lj elec e_short e_long e_corr 1/d\n");
      }
    }
    if(fph != NULL)  energy.iprint(istep,fph);
  }
}

void printSystem(int istep,SIMVARS &simvars,COORDS &coords)
{
  simvars.print();
  coords.print();
}
