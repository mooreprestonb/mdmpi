
#include "mdmol.h++"

static double *xm;

double getke(NEIGHBOR &ngbr)
{
  int dof = ngbr.getnmp()*DIM;  
  double *vm = xm+dof;
  double keng = .5*dot2(dof,vm);

#ifdef PARA
  double rkeng=keng;
  MPI_Reduce(&rkeng,&keng,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif

  return keng;
}

void average(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr,ENERGY &energy)
{
  double ke,pe,sr,lr,lj,elsr,ellr,elc;
  ke = getke(ngbr);
  pe = getpot(&sr,&lr,&lj,&elsr,&ellr,&elc);
  double d[DIM];
  cellDipole(coords.natoms(),coords.x(),coords.qch(),d);
  // double Q = cellQuadrupole(coords.natoms(),coords.x(),coords.qch());
  // printf("\nCell dipole = %g %g %g, d^2=%g, Q=%g\n",d[0],d[1],d[2],
  // (d[0]*d[0]+d[1]*d[1]+d[2]*d[2])*2.*M_PI/3.,2.*M_PI*Q/3.);
  energy.addAvg(ke,pe,sr,lr,lj,elsr,ellr,elc);
}

void startvv(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr)
{
  xm = ngbr.initMypart(coords);
  force(coords,ngbr);
#ifdef PARA
  int nmp = ngbr.getnmp();
  int dof = nmp*DIM;  
  getMyvals(nmp,coords.x(),xm,ngbr.getipmap());
  getMyvals(nmp,coords.v(),xm+dof,ngbr.getipmap());
  getMyvals(nmp,coords.a(),xm+dof*2,ngbr.getipmap());
#endif
  //int nmp = ngbr.getnmp();
  //double *am = xm+2*nmp*DIM;
  //for(int i=0;i<nmp;i++){ am[i] = 0.;}
  
}

void stepvv(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr)
{
  int i;  
  double dt = simvars.dt();
  double dt2 = .5*dt;
  int nmp = ngbr.getnmp();
  int dof = nmp*DIM;  
  double *vm = xm+dof;
  double *am = vm+dof;
  
  for(i=0;i<dof;i++) vm[i] += dt2*am[i];
  for(i=0;i<dof;i++) xm[i] += dt*vm[i];

  ngbr.update(coords);
  force(coords,ngbr);
#ifdef PARA
  getMyvals(nmp,coords.a(),am,ngbr.getipmap());
  // MPI_Barrier(MPI_COMM_WORLD); // lets wait for every proc to get here
#endif
  for(i=0;i<dof;i++) vm[i] += dt2*am[i];
}

void printcoords(void)
{
#ifdef VERBOSE
  for(i=0;i<nmp;i++) printf("%d x = %g %g %g, v= %g %g %g a= %g %g %g\n",i,
			    xm[i*DIM],xm[i*DIM+1],xm[i*DIM+2],
			    vm[i*DIM],vm[i*DIM+1],vm[i*DIM+2],
			    am[i*DIM],am[i*DIM+1],am[i*DIM+2]);
#endif
}
