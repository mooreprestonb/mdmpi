/*! \file neighbor.h++
  file contains the definition of the neighbor class
*/

#ifndef _NEIGHBOR_H
#define _NEIGHBOR_H

#include "coords.h++"

//! header for the neighbor class
class NEIGHBOR
{
private:
  int nintr,nix,nmol,nnow,nmp,nmpm,ipair,iperd,size,rank,ncell;
  int memvl,nvl,nexcl;
  int *ip,*ipvl,*hlist,*llist,*nlist,**exclij;
  int *ipmap,*disp,*counts,*icpmap;
  double *xi,*xdis,*recbuf,*xmpart;
  double skin,rcut;
  void checkSkin(COORDS &coords);
  void getVlist(double x[]);
  int nolist(void);
  // int nolist(double[]);
  int vl(void);
  //int vl(double x[]);
  int interLink(void);
  // int interLink(double x[]);
  int getLink(int ic,int ip[]);
  void getLnklist(double x[],int npnow);
  void updateLink(COORDS &coords);
  void recvcoords(double x[],double v[],int proc[]);
  void sendcoords(int np,int ipm[],int jcell,double x[],double v[]);
  void servReq(double x[]);
  int recvReq(int ipend[],double []);
  void recvcdsm(double x[],double v[],int proc[]);
  // void bcastIpmap(void);
  int packxi(int ,int ,int [],int [],double[],double[],int[],double[]);
  int packxi(int ,int ,int [],int [],int[]);
  int exclude(int,int);

public:
  NEIGHBOR(void); //!< constructor
  ~NEIGHBOR(void); //!< destrictor
  //! initialization with simulation vars and coordinates
  void init(SIMVARS &simvars,COORDS &coords);
  //! checks to see on needs to update the neighborlist 
  void update(COORDS &coords);
  //! get the interaction lists
  int interact(void);
  // int interact(double[]);
  //! sets the interaction list to -1 (ie restart the list)
  void setintr(void){nintr = -1;};
  //! returns the number of interactions that are in ip at the current state
  int getintrnow(void){return nnow;};
  //! returns point to interaction list
  int * getip(void){return ip;};
  //! returns point to interaction map
  int * getipmap(void){return ipmap;};
  //! returns pointer to temp buffer
  double * getxi(void){return xi;};
  //! returns pointer to receive buffer
  double * getrecbuf(void){return recbuf;};
  //! returns number of atoms in on this processor
  int getnmp(void){return nmp;};
  //! returns the type of pair list (ie the neighbor list type)
  int getipair(void){return ipair;};
  //! maps the coordinates onto this processors memory)
  double * initMypart(COORDS &coords);
  //! gathers particles onto centrial processor (ie for writing etc)
  void getcm(double [],double[]);
  //! add i and j to the exclucion list
  void addExclude(int i,int j);
  //! prints the exclusion list
  void printExclude(void);
};

#endif
