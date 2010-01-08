
#ifndef _GOR_IO_H
#define _GOR_IO_H
FILE *readHeader(char *name,int &natoms,int &nconf,double &dt);
int getConfig(FILE *fp,int natoms,double *x,double &cell);
#endif
