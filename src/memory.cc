
#include <stdlib.h>
#include <stdio.h>

//! routine to allocate memory and check that it was successful
void *cmalloc(size_t bytes)
{
  void *mem = malloc(bytes);
  if(mem==NULL){
    fprintf(stderr,"ERROR: can't allocate %d bytes (exiting)\n",(int)bytes);
    exit(1);
  }
  return mem;
}

void freec(void *mem){free(mem);}
