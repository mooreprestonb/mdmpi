/*! 
  \file coords.cc 
  \brief coords class functions 
  
  these function are part of the class external functions
*/

#include "coords.h++"
#include <sstream>
#include "xml.h++"
#include "parsexmll.h++"
#include "periodic.h++"

using namespace std;

/*! 
  \fn COORDS::COORDS(void)
  constuctor which sets
  _natoms to 0 
  and set array pointers to NULL
*/
COORDS :: COORDS(void)
{
  _natoms = _ntypes = 0;
  _itype = (int *) NULL;
  _icells = (int *) NULL;
  _x = _v = _a = (double *)NULL;
  _qch = _mass = (double *)NULL;
  for(int i=0;i<DIM;++i) _box[i] = -1.;
}

/*! 
  \fn COORDS::~COORDS(void)
  deconstuctor which frees arrays
*/
COORDS :: ~COORDS(void)
{
  if(_itype != (int *)NULL){free(_itype);}
  if(_x != (double *)NULL){free(_x);}
  if(_qch != (double *)NULL){free(_qch);}
}

#ifdef OLD
/*! 
  \fn oldread old style read
  \param name is the name of the input file
  \param natoms is the number of atoms
  \param x is the array of positions
  \param v is the array of velocities
  \param qch is the array of charges
  \param mass is the array of masses
*/
static int oldread(const string & name,int natoms,double *x,double *v,
		   double *qch,double *mass)
{
  const int MAXLINE = 256;
  int i,j,k,istep;
  char line[MAXLINE];
  FILE *fp;

  if((fp = fopen(name.c_str(),"r"))==NULL){
    fprintf(stderr,"ERROR: can't open coordfile \"%s\"\n",name.c_str());
    exit(1);
  }
  // read headder
  if(fgets(line,MAXLINE,fp)==NULL) {
    fprintf(stderr,"ERROR: empty file %s?\n",name.c_str());
    exit(1);
  } else {
    if(sscanf(line,"%*s %d %d",&i,&istep)!= 2){
      fprintf(stderr,"ERROR: parsing header of \"%s\"\n",name.c_str());
      exit(1);
    }
    if(natoms != i) {
      fprintf(stderr,"ERROR: # atoms do not match %d != %d\n",natoms,i);
      exit(1);
    }
  }
  
  for(i=0;i<natoms;++i){
    j = i*DIM;
    for(k=0;k<DIM;++k) 
      if(fscanf(fp,"%lg",&(x[j+k]))!=1) {
	fprintf(stderr,"ERROR: while reading in pos %d (dim %d out of %d)\n",
		i,k,DIM);
	exit(1);
      }
    for(k=0;k<DIM;++k) 
      if(fscanf(fp,"%lg",&(v[j+k]))!=1) {
	fprintf(stderr,"ERROR: while reading in vel %d (dim %d out of %d)\n",
		i,k,DIM);
	exit(1);
      }
    if(fscanf(fp,"%lg",&(qch[i])) != 1){
      fprintf(stderr,"ERROR: while reading in qch at %d\n",i);
      exit(1);
    }
    mass[i] = 1.;
//     if(fscanf(fp,"%lg",&(_mass[i])) != 1){
//       fprintf(stderr,"ERROR: while reading in qch at %d\n",i);
//       exit(1);
//     }
  }
  fclose(fp);
  return istep;
}
#endif

/*! 
  \brief function to read in the coordinates that are formated in 
  an xml type of manner

  \param name the name of the file to be parsed
  \param natoms the number of atoms
  \param *x a double array [DIM*natoms] that will contain the particles
  \param *v a double array [DIM*natoms] that will contain the velocities
  \param *qc a double array [natoms] that will contain the charges
  \param *mass a double array [natoms] that will contain the masses
  \param *itype a int array [natoms] that will contain the types of atoms
*/
static int readcml(const std::string & name,const int natoms,
		   double *x,double *v,double *qch,double *mass,int *itype)
{
  char line[256];
  sprintf(line,"Reading in coords from file %s\n",name.c_str());
  writeLog(line);
  
  list<XML> lxml = parseXMLfile(name);
  XML txml = lxml.front();
  string stmp = "coords";
  XML txml2 = txml.findtag(stmp);
  if(txml == txml2) {
    cerr << "ERROR: coords tag not found in file \""<<name<<"\"\n";
    exit(1);
  }

  XML txml3;
  txml3 = txml2.findtag("natoms");
  if(txml3 == txml2) {
    cerr << "ERROR: natoms tag not found in file \""<<name<<"\"\n";
    exit(1);
  }
  int npart = atoi(txml3.data().front().c_str());
  if(npart != natoms) {
    cerr << "ERROR: natoms don't match "<<npart<<" != "<<natoms<<endl;
    exit(1);
  }
  int i,j;
  list<XML>::iterator it,end;
  list<XML> cxml = txml2.xmllist();
  end = cxml.end();
  i=0;
  for(it=cxml.begin();it != end;++it){
    if((*it).tag() == "xyz") {
      j = i*DIM;
      istringstream is((*it).data().front());
      is >> x[j] >> x[j+1] >> x[j+2];
      ++i;
    }
  }
  i=0;
  for(it=cxml.begin();it != end;++it){
    if((*it).tag() == "vxyz") {
      j = i*DIM;
      istringstream is((*it).data().front());
      is >> v[j] >> v[j+1] >> v[j+2];
      ++i;
    }
  }
  i=0;
  for(it=cxml.begin();it != end;++it){
    if((*it).tag() == "charge") {
      istringstream is((*it).data().front());
      is >> qch[i];
      ++i;
    }
  }
  i=0;
  for(it=cxml.begin();it != end;++it){
    if((*it).tag() == "mass") {
      istringstream is((*it).data().front());
      is >> mass[i];
      ++i;
    }
  }
  i=0;
  for(it=cxml.begin();it != end;++it){
    if((*it).tag() == "type") {
      istringstream is((*it).data().front());
      is >> itype[i];
      ++i;
    }
  }
  txml3 = txml2.findtag("step");
  if(txml3 == txml2) exit(1);
  int istep = atoi(txml3.data().front().c_str());
  return istep;
}


/*!
  Set the mass to 1 if not set
 */
void setmass(int natoms,double mass[])
{ 
  for(int i=0;i<natoms;++i) {
    if(mass[i]==0) {
      mass[i] = 1.;
    }
  }
}

/*! 
  \fn int COORDS::read(int n,const std::string &name)
  \brief reads coordinates from file name \a name the format is

  \param n number of atoms in the coordinate file
  this parameter is for consistancy...
  \param name a const string that is the name of the coordinate file
<p>
*/
int COORDS :: read(int natoms,const std::string & name)
{
  int istep;
  _natoms = natoms;
  if(_x != NULL) free(_x);
  int dof = natoms*DIM;
  _x = (double *)calloc(3*dof,sizeof(double));
  _v = _x + dof;
  _a = _x + 2*dof;
  // for(i=0;i<3*dof;++i) x[i] = 0.;
  _qch = (double *)calloc(2*natoms,sizeof(double));
  _mass = _qch +natoms;

  _ntypes = 1;
  _itype = (int *)calloc(natoms,sizeof(int));
  // read in file 
#ifdef OLD
  istep = oldread(name,natoms,_x,_v,_qch,_mass);
#else
  istep = readcml(name,natoms,_x,_v,_qch,_mass,_itype);
#endif
  setmass(natoms,_mass);
  return istep;
}

void COORDS::reset(void)
{
  resetcoords(_natoms,_x,_v);
}

/*! 
  \fn void COORDS::print(FILE *fp=stdout)
  prints the coords to a \a fp or stdout if non is given
  \param fp file pointer
*/
void COORDS::print(FILE *fp)
{
  int i,j,k;
  fprintf(fp,"<coords> <natoms> %d </natoms>\n",_natoms);
  for(i=0;i<_natoms;++i){
    j = DIM*i;
    fprintf(fp,"<xyz id=\"%d\"> ",i);
    for(k=0;k<DIM;++k) fprintf(fp,"%g ",_x[j+k]); 
    fprintf(fp,"</xyz>\n<vxyz id=\"%d\"> ",i);
    for(k=0;k<DIM;++k) fprintf(fp,"%g ",_v[j+k]); 
    fprintf(fp,"</vxyz>\n<axyz id=\"%d\">",i);
    for(k=0;k<DIM;++k) fprintf(fp,"%g ",_a[j+k]); 
    fprintf(fp,"</axyz>\n");
    fprintf(fp,"<charge id=\"%d\"> %g </charge> ",i,_qch[i]);
    fprintf(fp,"<mass id=\"%d\"> %g </mass>\n",i,_mass[i]);
  }
  fprintf(fp,"</coords>");
}
