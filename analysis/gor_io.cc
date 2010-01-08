

#include <stdio.h>
#include <stdlib.h>
#include "gor_io.h++"

const int DIM=3;
const int MAXLINELEN=256;

FILE *readHeader(char *name,int &natoms,int &nconf,double &dt)
{
  FILE *fp;
  /* open configuration file  (check to see that it opened correctly) */
  if((fp = fopen(name,"r")) == NULL){
    fprintf(stderr,"ERROR: can't open %s (conf-file)\n",name);
    exit(1);
  }

  /* read in header info */
  char line[MAXLINELEN];
  if(fgets(line,MAXLINELEN,fp) == NULL){
    fprintf(stderr,"ERROR: while reading in header\n");
    exit(1);
  }
  
  sscanf(line,"%*s %d %d %lg\n",&natoms,&nconf,&dt);
  return fp;
}
/*----------------------------------------------------------*/
int getConfig(FILE *fp,int natoms,double *x,double &cell)
{
  int i;
  char line[MAXLINELEN];
  
  if(fgets(line,MAXLINELEN,fp) == NULL){
    /* you hit the EOF so return to kick out of while loop */
    return 0;
  }

  /* parse up this line */
  if(sscanf(line,"%lg %lg %lg",&x[0],&x[1],&x[2]) != 3){
    fprintf(stderr,"ERROR: something is wrong with configuration file\n");
    exit(1);
  }
  
  /* loop over rest of atoms and make sure we read in what we should */
  for(i=1;i<natoms;i++){
    if(fgets(line,MAXLINELEN,fp) == NULL){
      fprintf(stderr,"ERROR: something is wrong with configuration file\n");
      exit(1);
    }
    if(sscanf(line,"%lg %lg %lg",&x[i*DIM],&x[i*DIM+1],&x[i*DIM+2]) != 3){
      fprintf(stderr,"ERROR: something is wrong with configuration file\n");
      exit(1);
    }
  }
  /* read in the cell information */
  if(fgets(line,MAXLINELEN,fp) == NULL){
    fprintf(stderr,"ERROR: something is wrong with configuration file\n");
    exit(1);
  } /* should put in error checking here so that we make sure it is diagonal */
  if(sscanf(line,"%lg %*g %*g %*g %*g %*g %*g %*g %*g",&cell) != 1){
    fprintf(stderr,"ERROR: something is wrong with configuration file\n");
    exit(1);
  }
  return 1;  /* still more to read */
}  
