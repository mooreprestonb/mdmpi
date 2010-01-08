
#ifndef _MD_IO_H
#define _MD_IO_H
void writeELog(int istep,SIMVARS &simvars,ENERGY &energy);
void saveRestart(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr,ENERGY &);
void storeConf(SIMVARS &simvars,COORDS &coords,NEIGHBOR &ngbr);
void storeHam(const char *filen,int istep,ENERGY &energy,int rank);
void printSystem(int,SIMVARS &simvars,COORDS &coords);
void writeOut(const char *s);
void writeErr(const char *s);
void setStderr(const char *file);
void setStdout(const char *file);
FILE * setLogfile(const char *file);
void writeLog(const char *s);
#endif
