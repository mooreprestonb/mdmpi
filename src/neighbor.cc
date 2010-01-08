
#include "mdmol.h++"

const int NIX = 501;

/*-------------------------------------------------------------*/
NEIGHBOR :: NEIGHBOR(void){
  nix=NIX;nexcl=0;
  ipair=nmol=ncell,nintr=nnow=nmpm=nmp=nvl=memvl=size=rank=iperd=0;
  ip=hlist=llist=nlist=ipmap=ipvl=disp=counts=icpmap=NULL;
  xi=xdis=recbuf=xmpart=NULL;
  exclij=NULL;
}

/*-------------------------------------------------------------*/
NEIGHBOR :: ~NEIGHBOR(void){
  if(ip != NULL) delete[] ip;
  if(ipvl != NULL) delete[] ipvl;
  if(ipmap != NULL) delete[] ipmap;
  if(hlist != NULL) delete[] hlist;
  if(nlist != NULL) delete[] nlist;
  if(llist != NULL) delete[] llist;
  if(disp != NULL) delete[] disp;
  if(counts != NULL) delete[] counts;
  if(xi != NULL) delete[] xi;
  if(xdis != NULL) delete[] xdis;
  if(recbuf != NULL) delete[] recbuf;
  if(icpmap != NULL) delete[] icpmap;
  if(exclij!=NULL){
    for(int i=0;i<nmol;i++) free((int *)exclij[i]);
    free(exclij);
  }
}

void NEIGHBOR :: printExclude(void)
{
  char line[MAXLINE];
  int i,k;
  sprintf(line,"# Exclusions = %d\n",nexcl); writeErr(line);
  for(i=0;i<nmol;i++){
    int *eij = exclij[i];
    if(eij[0] != -1) {
      sprintf(line,"Exclude %d from ",i);
      k = 0;
      while(eij[k] != -1) sprintf(line,"%s%d ",line,eij[k++]);
      writeErr(line);
    }
  }
}

void NEIGHBOR :: addExclude(int i,int j)
{
  int k;
  if(j<i){k=i;i=j;j=k;}
  int *eij = exclij[i];
  k=0;
  while(eij[k] != -1) {if(eij[k++]==j) return;}
  nexcl++;
  exclij[i] = (int *)realloc(exclij[i],(k+2)*sizeof(int));
  exclij[i][k] = j;
  exclij[i][k+1] = -1;
  if(iperd==3) addEcorr(i,j);  
}

int NEIGHBOR :: exclude(int i,int j)
{
  int k;
  if(i==j) return 1;
  if(j<i){k=i;i=j;j=k;}
  int *eij = exclij[i];
  k=-1;
  while(eij[++k] != -1) if(eij[k]==j) return 1;
  return 0;
} 

/*-------------------------------------------------------------*/
void NEIGHBOR :: init(SIMVARS &simvars,COORDS &coords)
{
  char line[MAXLINE];
  size = simvars.size();
  rank = simvars.rank();
  nmol = simvars.nmol();
  iperd = simvars.iperd();
  ip = new int[2*nix];
  xi = new double[2*nix*DIM];
  ncell = simvars.ncell();
  ipair = simvars.ipair();
  rcut  = simvars.rcut();
  skin  = simvars.skin();

  exclij = (int **)malloc(nmol*sizeof(int *));
  int i;
  for(i=0;i<nmol;i++){
    exclij[i] = (int *)malloc(sizeof(int));
    exclij[i][0] = -1;
  }
  switch(ipair){ /* initial neigbor lists */
  case 0: case 1:
    if(rank==0) { 
      if(ipair==0) sprintf(line,"N^2 algorithm being used.");
      else sprintf(line,"N^2 algorithm2 being used.");
      writeOut(line);
    }
    break;  
  case 2:   case 3:
    if(rank==0) {
      if (ipair == 2) sprintf(line,"Verlist being used.");
      else sprintf(line,"Verlist2 being used.");
    }
    writeOut(line);
    llist = new int[nmol];
    checkSkin(coords);
    break;
  case 4: case 5: 
    if(rank==0) {
      if(ncell<3){
	sprintf(line,"Link list being used but ncell=%d<3",ncell);
	writeErr(line);
	exit(1);
      }
      if(ipair == 4) sprintf(line,"Link list being used ncell=%d ~ %gA.",
			     ncell,simvars.xbox()/ncell);
      else sprintf(line,"Link2 list being used ncell=%d ~ %gA.",
		  ncell,simvars.xbox()/ncell);
      writeOut(line);
    }
    hlist = new int[ncell*ncell*ncell];
    nlist = new int[ncell*ncell*ncell];
    icpmap = new int[ncell*ncell*ncell];
    int i;
    for(i=0;i<size;i++){
      int ib,ie,j;
      decomp1d(ncell*ncell*ncell,size,i,&ib,&ie); 
      for(j=ib;j<ie;j++) icpmap[j] = i;
    }
    llist = new int[nmol];
    break;
  default:
    sprintf(line,"ERROR: in get_neigbors ipair=%d\n",ipair); 
    writeErr(line);
    exit(1);
  }
}

/*-------------------------------------------------------------*/
void NEIGHBOR :: update(COORDS &coords)
{
  int nmpt;
  
  nintr=-1; // reset to get all the list
  nmpt=nmp; //nmp needs to be part of NEIGHBOR?
  
  switch(ipair){ /* initial neigbor lists */
  case 0:  case 1: 
    break;
  case 2:  case 3: 
    checkSkin(coords);    
    break;
  case 4:  case 5: updateLink(coords); break;
  default:
    fprintf(stderr,"ERROR: in ud_neigbors ipair=%d\n",ipair); exit(1);
  }
}

/*--------------------------------------------------------*/
int NEIGHBOR :: interact(void)
{
  int done=0;

  switch(ipair){
  case 0: case 1: done = nolist(); break;
  case 2: case 3: done = vl();  break;
  case 4: case 5: done = interLink(); break;
  default:fprintf(stderr,"ERROR: in get_neigbors ipair=%d\n",ipair); exit(1);
  }
  return done;
}

/*--------------------------------------------------------*/
double * NEIGHBOR :: getrecbuf(int n){
  if (recbuf != NULL){
    delete recbuf;
  }
  recbuf = new double[n];
  return recbuf;
};

/*--------------------------------------------------------*/
