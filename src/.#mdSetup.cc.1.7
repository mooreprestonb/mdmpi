/*! \file mdSetup.cc
  file that sets up the simulation variables
*/

#include "mdmol.h++"

//! read input and assign variables...
void setup(FLAGS &flags,SIMVARS &simvars,COORDS &coords,
	   NEIGHBOR &ngbr,ENERGY &energy)
{
  char name[MAXLINE],line[MAXLINE];
  int rank = 0;
  int size = 1;
  gethostname(name,MAXLINE);

  simvars.hostname(name);
  simvars.pid((int)getpid());

  simvars.inname(flags.infile());
  
#ifdef PARA
  int argc = flags.argc();
  char **argv = flags.argv();
  MPI_Init(&argc, &argv);  /* initialize MPI */
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  fprintf(stdout,"Running MPI Version %s on %s (id=%d) of rank %d and size %d with input %s\n",
          argv[0],name,simvars.pid(),rank,size,flags.infile());
#ifdef SLEEP
  printf("Sleeping for 30 secs\n");
  MPI_Barrier(MPI_COMM_WORLD);
  sleep(30); // sleep so that we can debug
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#else
  fprintf(stdout,"Running %s on %s (id=%d)\n",flags.command(),name,simvars.pid());
#endif
  simvars.rank(rank); 
  simvars.size(size);

  if(flags.q()){
    char hostname[256];
    gethostname(hostname,256);
    fprintf(stdout,"Running %s on %s (id=%d), reading input %s\n",
	    flags.command(),hostname,simvars.pid(),flags.infile());
    // fprintf(stdout,"infile = \"%s\"\n",_infile);
  }



  /* setup */
  simvars.readinit();

  int istep = coords.read(simvars.nmol(),simvars.coordname());
  if(simvars.calctype()==1){
    energy.read(simvars.nmol(),simvars.coordname());
    simvars.istep(istep);
  } else {
    simvars.istep(0); // don't need this but lets be safe
  }

  ewaldSetup(simvars.nmol(),coords.qch(),simvars.kmax(),simvars.alpha(),
	     simvars.rcut(),simvars.esurf(),simvars.iperd(),simvars.cutoff());

  { // LJ setup
    int ntypes = 1;
    double sig[1],eps[1];
    sig[0] = eps[0] = 1.;
    LJsetup(ntypes,eps,sig,simvars.rcut(),simvars.xbox());
  }

  setIbox(simvars.iperd(),simvars.xbox()); // set periodic conditions
  ngbr.init(simvars,coords); // initiallize neighbor lists

  setupFreeze(simvars.nfreeze()); // setup freeze list
  for(int i=0;i<simvars.nfreeze()-1;++i){
    for(int j=i+1;j<simvars.nfreeze();++j) ngbr.addExclude(i,j);
  }

  if(simvars.rank()==0) {
    sprintf(line,"Getting initial energy and force");
    writeOut(line);
  }

  startvv(simvars,coords,ngbr); // initialize integrator

  energy.dof(simvars.nmol()*DIM);

  if(simvars.rank()==0) {
    sprintf(line,"Starting MD ntime = %d with %d dof",
	    simvars.nstep(),simvars.dof());
    writeOut(line);
  }
}
