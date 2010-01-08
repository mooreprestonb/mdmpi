
enum KEYWORD {
  NATOMS,   NSTEPS,   TIMESTEP, BOX,     COORDFILE,  
  NEIGHBOR, PERIODIC, NCELL,    RCUT,    SKIN, 
  KMAX,     ALPHA,    ESURF,    WRITE,   NFREEZE,
  RESTART,  TRAJ,     ENERGY,   LOG,     ERR,
  CALC,     TMZERO,   CUTOFF,   INITRUN, CONTRUN
};

char keyword[25][50] = {
  "natoms",  "nsteps",  "timestep",  "box",  "coordsfile",
  "neighbor",  "periodicity",  "ncell",  "rcut",  "skin",
  "kmax",  "alpha",  "esurf",  "writefreq",  "nfreeze",
  "restartfile",  "trajectoryfile",  "energyfile",  "logfile",  "errorfile",
  "calculation",  "zerotm",  "cutoff",  "initial",  "continue"
};
