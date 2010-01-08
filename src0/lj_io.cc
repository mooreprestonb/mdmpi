
#include "mdmol.h++"


void save_rest(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr)
{
  char *fn;
  int i,j,k,nmol;
  FILE *fp;

  double *x = coords.getx();
  double *v = coords.getv();    
  ngbr.getcm(x,v);

  if(simvars.rank==0){
    nmol = simvars.nmol;
    i = strlen(simvars.getrestname());
    fn = (char *)malloc((i+10)*sizeof(char));
    sprintf(fn,"%s.%d",simvars.getrestname(),simvars.ntime);
    printf("Writing restart file %s\n",fn);
    if((fp=fopen(fn,"w"))==NULL){
      fprintf(stderr,"ERROR: can't open %s\n",fn);
      exit(1);
    }
    free(fn);
#ifdef ZEROTM
    {
      int dof=nmol*DIM;
      double *tv = (double *)malloc(DIM*sizeof(double));
      double rn = 1./(double)nmol;
      double *v = coords.getv();
      for(k=0;k<DIM;k++) tv[k]=0.;
      for(i=0;i<dof;i+=DIM) for(k=0;k<DIM;k++) tv[k] += v[i+k];
      for(k=0;k<DIM;k++) tv[k] *= rn;
      for(i=0;i<dof;i+=DIM) for(k=0;k<DIM;k++) v[i+k] -= tv[k];
      free(tv);
    }
#endif

#ifdef RESETPOS
    periodic(nmol,x);
#endif
    fprintf(fp,"# %d %d %g %g %d %d %d ",nmol,simvars.ntime,simvars.dt,
	    simvars.xbox,simvars.ipot,DIM,simvars.iperd);
    fprintf(fp,"%s %s %s\n",simvars.getrestname(),simvars.gethamname(),
	    simvars.getconfname());
    
    for(i=0;i<nmol;i++){
      j = i*DIM;
      for(k=0;k<DIM;k++) fprintf(fp,"%g ",x[j+k]);
      fprintf(fp,"\t");
      for(k=0;k<DIM;k++) fprintf(fp,"%g ",v[j+k]);
      fprintf(fp,"\n");
    }
#ifdef WRITE_BOX
    fprintf(fp,"%g 0 0  0 %g 0  0 0 %g\n",
	    simvars.xbox,simvars.xbox,simvars.xbox);
#endif
    fclose(fp);
  }
}

void store_conf(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr)
{
  static int istart=0;
  static FILE *fp;
  int nmol,i,j,k;

  double *x = coords.getx();
  double *v = coords.getv();
  ngbr.getcm(x,v);

  if(simvars.rank==0){
    nmol = simvars.nmol;
  
    if(istart==0){
      /* printf("Enter config name:"); scanf("%s",filen); */
      if((fp = fopen(simvars.getconfname(),"w"))==NULL){
	fprintf(stderr,"ERROR:can't open %s\n",simvars.getconfname());
	exit(1);
      }
      istart=1;
      fprintf(fp,"# %d %g %g %d %d %d\n",nmol,simvars.dt,
	      simvars.xbox,NCONF,simvars.ipot,simvars.iperd);
    }

    for(i=0;i<nmol;i++){
      j = DIM*i;
      for(k=0;k<DIM;k++) fprintf(fp,"%g ",x[j+k]);
      for(k=DIM;k<WDIM;k++) fprintf(fp,"0 ");
      fprintf(fp,"\n");
    }
#ifdef WRITE_BOX
    fprintf(fp,"%g 0 0  0 %g 0  0 0 %g\n",
	    simvars.xbox,simvars.xbox,simvars.xbox);
#endif
  }
}

void store_ham(char *filen,int istep,ENERGY &energy,int rank)
{
  static int istart=0;
  static FILE *fp;
  double rd = 1./energy.dof();
  
  if(rank==0){
    if(istart==0){
      istart=1;
      if((fp = fopen(filen,"w"))==NULL){
	fprintf(stderr,"ERROR: can't open %s\n",filen);
	exit(1);
      }
      fprintf(fp,"# istep ham ke pe dof=%d\n",energy.dof());
    }
    fprintf(fp,"%d %g %g %g\n",istep,energy.ham()*rd,
	    energy.ke()*rd,energy.pe()*rd);
  }
}

void print_system(SIMVARS &simvars,COORDS &coords)
{
  int i,j,k;

  double *x = coords.getx();
  double *v = coords.getv();
  double *a = coords.geta();
  
  printf("rank=%d nmol=%d size=%d\n",simvars.rank,simvars.nmol,simvars.size);
  for(i=0;i< simvars.nmol;i++){
    j = DIM*i;
    printf("%d %d ",simvars.rank,i);
    for(k=0;k<DIM;k++) printf("%6.3g ",x[j+k]);
    for(k=0;k<DIM;k++) printf("%6.3g ",v[j+k]);
    for(k=0;k<DIM;k++) printf("%6.3g ",a[j+k]);
    printf("\n");
    fflush(stdout);
  }
}
