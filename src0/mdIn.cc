/*! \file mdIn.cc
  file to parse and readin the simulation input params
*/

#include <stdio.h>
#include <iostream>
#include <list>
#include <algorithm>
#include <string>
#include <sstream>

#include "mdmol.h++"
#include "simkeys.h++"
#include "xml.h++"

#include "parsexmll.h++"

using namespace std;

string keyword[27] = {
  "natoms",     "nsteps",        "timestep",  "box",      "coordsfile",
  "neighbor",   "periodicity",   "ncell",     "rcut",     "skin",
  "kmax",       "alpha",         "esurf",     "writefreq","nfreeze",
  "restartfile","trajectoryfile","energyfile","logfile",  "errorfile",
  "calculation","zerotm",        "cutoff",    "initial",  "continue",
  "systemfile", "paramfile"
};

const int NKEYWORD=27;

enum KEYWORD {
  NATOMS,   NSTEPS,   TIMESTEP, BOX,     COORDFILE,  
  NEIGHBOR, PERIODIC, NCELL,    RCUT,    SKIN, 
  KMAX,     ALPHA,    ESURF,    WRITE,   NFREEZE,
  RESTART,  TRAJ,     ENERGY,   LOG,     ERR,
  CALC,     TMZERO,   CUTOFF,   INITRUN, CONTRUN,
  SYSTEM,   PARAM
};

string systemtype[] = {"system","molecule","name","index","nmol"};
string forcewords[] = {"forceparm","inter","atom1","atom2",
		  "pottype","sigma","epsilon"};

bool atob(const string & value)
{
  if(value[0] == 't' || value[0] == 'T') return true;
  if(value == "on") return true;
  if(value == "ON") return true;
  return false;
}

string nospacestring(const string &src)
{
  istringstream is(src);
  string out;
  is >> out;
  return out;
}

/*! sets the default keys and initial values
  
*/
void setkeys(list<SIMKEY> &keys)
{
  SIMKEY key;

#ifdef NOSTRING
  key.setdef(1,int(SIMKEY::STRING),keyword[SYSTEM].c_str(),NULL); keys.push_back(key);
  key.setdef(1,int(SIMKEY::INT),keyword[NATOMS].c_str(),NULL);   keys.push_back(key);
  key.setdef(1,int(SIMKEY::INT),keyword[NSTEPS].c_str(),NULL);   keys.push_back(key);
  key.setdef(1,int(SIMKEY::DOUBLE),keyword[TIMESTEP].c_str(),NULL);  keys.push_back(key);
  key.setdef(1,int(SIMKEY::DOUBLE),keyword[BOX].c_str(),NULL);   keys.push_back(key);
  key.setdef(1,int(SIMKEY::STRING),keyword[COORDFILE].c_str(),NULL);  keys.push_back(key);
  key.setdef(1,int(SIMKEY::STRING),keyword[PARAM].c_str(),NULL);   keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[NEIGHBOR].c_str(),"2");   keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[PERIODIC].c_str(),"3");   keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[NCELL].c_str(),"3");   keys.push_back(key);
  key.setdef(0,int(SIMKEY::DOUBLE),keyword[RCUT].c_str(),"2.5"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::DOUBLE),keyword[SKIN].c_str(),"1.0"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[KMAX].c_str(),"10"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::DOUBLE),keyword[ALPHA].c_str(),"1.0"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::DOUBLE),keyword[ESURF].c_str(),"-1.0"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[WRITE].c_str(),"10"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[NFREEZE].c_str(),"0"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[RESTART].c_str(),NULL); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[TRAJ].c_str(),NULL); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[ENERGY].c_str(),NULL); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[LOG].c_str(),NULL); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[CALC].c_str(),keyword[INITRUN].c_str()); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[TMZERO].c_str(),"on"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[CUTOFF].c_str(),"off"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[ERR].c_str(),NULL); keys.push_back(key);
#else
  string null = "";
  key.setdef(1,int(SIMKEY::STRING),keyword[SYSTEM],null); keys.push_back(key);
  key.setdef(1,int(SIMKEY::INT),keyword[NATOMS],null);   keys.push_back(key);
  key.setdef(1,int(SIMKEY::INT),keyword[NSTEPS],null);   keys.push_back(key);
  key.setdef(1,int(SIMKEY::DOUBLE),keyword[TIMESTEP],null);  keys.push_back(key);
  key.setdef(1,int(SIMKEY::DOUBLE),keyword[BOX],null);   keys.push_back(key);
  key.setdef(1,int(SIMKEY::STRING),keyword[COORDFILE],null);  keys.push_back(key);
  key.setdef(1,int(SIMKEY::STRING),keyword[PARAM],null);   keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[NEIGHBOR],"2");   keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[PERIODIC],"3");   keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[NCELL],"3");   keys.push_back(key);
  key.setdef(0,int(SIMKEY::DOUBLE),keyword[RCUT],"2.5"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::DOUBLE),keyword[SKIN],"1.0"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[KMAX],"10"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::DOUBLE),keyword[ALPHA],"1.0"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::DOUBLE),keyword[ESURF],"-1.0"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[WRITE],"10"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::INT),keyword[NFREEZE],"0"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[RESTART],null); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[TRAJ],null); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[ENERGY],null); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[LOG],null); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[CALC],keyword[INITRUN]); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[TMZERO],"on"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[CUTOFF],"off"); keys.push_back(key);
  key.setdef(0,int(SIMKEY::STRING),keyword[ERR],null); keys.push_back(key);
#endif
}

void keymerge(list<SIMKEY> &keys,list<SIMKEY> &dkeys)
{
  int ierr=0;
  int inot=0;
  list<SIMKEY>::iterator iter1,iter2;
  
  // make sure we have the required key words 
  for(iter1 = dkeys.begin();iter1 != dkeys.end();++iter1){
    iter2=find(keys.begin(),keys.end(), (*iter1).key());
    if((*iter1).isReq() && iter2==keys.end()) {
#ifdef NOSTRING
      fprintf(stderr,"Required keyword \"%s\" not found\n",(*iter1).key());
#else
      fprintf(stderr,"Required keyword \"%s\" not found\n",(*iter1).key().c_str());
#endif
      ierr++;
    } else { // we match so merge key, then erase set ...
      if(iter2 != keys.end()) {
	(*iter1) += (*iter2);
	keys.erase(iter2);
      } 
    }
  } 
  // Keys that we don't recognize 
  if(keys.begin() != keys.end()){
    fprintf(stderr,"The following %d keywords were not recognized\n",inot);
    for(iter1 = keys.begin();iter1 != keys.end();++iter1){
      (*iter1).print(stderr);
      ierr++;
    }
  }
  
  for(iter1 = dkeys.begin();iter1 != dkeys.end();++iter1){
    if((*iter1).iset() > 1){
#ifdef NOSTRING
      fprintf(stderr,"keyword \"%s\" found %d time!\n",
	      (*iter1).key(),(*iter1).iset());
#else
      fprintf(stderr,"keyword \"%s\" found %d time!\n",
	      (*iter1).key().c_str(),(*iter1).iset());
#endif
      // ierr += (*iter1).iset()-1;
      ++ierr;
    }
  }
  if(ierr!=0) {
    fprintf(stderr,"There were %d errors in the input\n",ierr);
    fprintf(stderr,"We have the following keys (with their present values)\n");
    for(iter1=dkeys.begin();iter1!=dkeys.end();++iter1) cout<<*iter1<<endl;
    exit(ierr);
  }
}

void getkeys(list<XML> &lxml,list<SIMKEY> &keys)
{
  string w1,w2;
  SIMKEY key;
  XML txml,txml2;
  list<XML>::iterator it,end;
  txml = lxml.front();
  txml2 = txml.findtag("simulation");
  if(txml == txml2) {
    cerr << "key \"simulation\" not found\n";
    exit(1);
  }
  it = txml2.xmllist().begin();
  end = txml2.xmllist().end();
  for(;it != end;++it){
    string w1 = nospacestring((*it).tag());
    string w2 = nospacestring((*it).data().front());
#ifdef NOSTRING
    key.set(w1.c_str(),w2.c_str());
#else
    key.set(w1,w2);
#endif
    keys.push_back(key);
  }
}

int getsnline(FILE *fp,int n,char *name,SIMKEY &key)
{
  int done = 0;
  char line[MAXLINE];
  char w1[MAXLINE];
  char w2[MAXLINE];
  while(!done){
    if(fgets(line,MAXLINE,fp)==NULL) {
      //fprintf(stderr,"ERROR: at EOF on line %d of %s?\n",n,name);
      //exit(1);
      return 0;
    } else {
      if(sscanf(line,"%s %s",w1,w2)!= 2){
 	//fprintf(stderr,"ERROR: reading in line %d of \"%s\"\n",n,name);
 	//exit(1);
 	// return 0; 
      } else {
 	if(w1[0] != '#') {
 	  done = 1;
 	  break; // we could say done... but we break instead;
 	}
      }
    }
  }
  key.set(w1,w2);
  return 1;
}

void oldgetkeys(char *setname,list<SIMKEY> &keys)
{
  int nline=0;
  FILE *fp;
  SIMKEY key;
  if((fp = fopen(setname,"r"))==NULL){
    fprintf(stderr,"ERROR: can't open old input file \"%s\"\n",setname);
    exit(1);
  }
  // read in keys
  while(getsnline(fp,++nline,setname,key)) keys.push_back(key);
  fclose(fp);
}   

#define SHOWKEYS
void readInput(string &setname, list<SIMKEY> & dkeys)
{
  SIMKEY key;
  list<SIMKEY> keys;
  
  setkeys(dkeys); //set default list

#ifdef SHOWKEYS
  {
    list<SIMKEY>::iterator iter;
    for(iter = dkeys.begin();iter!=dkeys.end();++iter) cout << *iter <<endl;
  }
#endif

  // read in keys
#ifdef OLD
    oldgetkeys(setname,keys);
#else
    list<XML> lxml = parseXMLfile(setname);
    getkeys(lxml,keys);
#endif
  
#ifdef SHOWKEYS
  {
    list<SIMKEY>::iterator iter;
    for(iter = keys.begin();iter!=keys.end();++iter) cout << *iter <<endl;
  }
#endif

  keymerge(keys,dkeys);
}

void SIMVARS::readinit(void)
{
  if(_rank != 0) return;
  
  int i;
  char line[MAXLINE];
  list<SIMKEY> keys;
  list<SIMKEY>::iterator iter;

  readInput(_inname,keys);

  for(i=0;i<NKEYWORD;++i){
    iter = find(keys.begin(),keys.end(),keyword[i]);
    if(iter == keys.end()) {
      fprintf(stderr,"Keyword \"keyword[%d] = %s\"not found check dictionary!!!\n",
	      i,keyword[i].c_str());
      continue;
    }
    switch(i){
    case NATOMS: 
      _nmol = atoi((*iter).value().c_str());  
      break;
    case NSTEPS:   
      _nstep = atoi((*iter).value().c_str()); 
      break;
    case TIMESTEP: 
      _dt = atof((*iter).value().c_str()); 
      break;
    case BOX: 
      _xbox = atof((*iter).value().c_str()); 
      break;
    case COORDFILE: 
      coordname((*iter).value()); 
      break;
    case NEIGHBOR: 
      _ipair = atoi((*iter).value().c_str()); 
      break; 
    case PERIODIC: 
      _iperd = atoi((*iter).value().c_str()); 
      break;
    case NCELL: 
      _ncell = atoi((*iter).value().c_str()); 
      break;
    case RCUT: 
      _rcut = atof((*iter).value().c_str()); 
      break;
    case SKIN: 
      _skin = atof((*iter).value().c_str()); 
      break;
    case KMAX: 
      _kmax = atoi((*iter).value().c_str()); 
      break;
    case ALPHA: 
      _alpha = atof((*iter).value().c_str());
      break;
    case ESURF:
      _esurf = atof((*iter).value().c_str());
      break;    
    case WRITE: 
      _iwrite = atoi((*iter).value().c_str()); 
      break;
    case NFREEZE: 
      _nfreeze = atoi((*iter).value().c_str()); 
      break;
    case RESTART: 
      restname((*iter).value()); 
      break;
    case TRAJ: 
      trajname((*iter).value()); 
      break;
    case ENERGY: 
      hamname((*iter).value()); 
      break;
    case LOG:
      logname((*iter).value()); 
      break;
    case ERR: 
      errname((*iter).value()); 
      break;
    case CALC: 
      if((*iter).value()== keyword[INITRUN]) _calctype = 0;
      else if((*iter).value() == keyword[CONTRUN]) _calctype = 1;
      else {
	sprintf(line,"Error in input file \"%s\" calculation type = \"%s\"",
		_inname.c_str(),(*iter).value().c_str()); writeErr(line);
	sprintf(line,"Can be either \"%s\" or \"%s\"\n",
		keyword[INITRUN].c_str(),keyword[CONTRUN].c_str()); 
	writeErr(line);
	exit(1);
      }
      break;
    case TMZERO:  _zerotm = atob((*iter).value()); break;
    case CUTOFF:   _cutoff = atob((*iter).value());break;
    case SYSTEM: systemfile((*iter).value());       break;
    case PARAM: 
      paramfile((*iter).value()); 
      break;
    default:
      fprintf(stderr,"Error???\n");
      (*iter).print(stderr);
      break;
    }
    keys.erase(iter);
  }

  
#ifdef NOSTRING
  iter = find(keys.begin(),keys.end(),keyword[NEIGHBOR].c_str());
  _ipair = atoi((*iter).value()); 
  // iperd
  iter = find(keys.begin(),keys.end(),keyword[PERIODIC].c_str());
  _iperd = atoi((*iter).value()); 
  // ncell 
  iter = find(keys.begin(),keys.end(),keyword[NCELL].c_str());
  _ncell = atoi((*iter).value()); 
  // rcut 
  iter = find(keys.begin(),keys.end(),keyword[RCUT].c_str());
  _rcut = atof((*iter).value()); 
  // skin
  iter = find(keys.begin(),keys.end(),keyword[SKIN].c_str());
  _skin = atof((*iter).value()); 
  
  // kmax 
  iter = find(keys.begin(),keys.end(),keyword[KMAX].c_str());
  _kmax = atoi((*iter).value()); 
  // alpha 
  iter = find(keys.begin(),keys.end(),keyword[ALPHA].c_str());
  _alpha = atof((*iter).value());
  // esurf 
  iter = find(keys.begin(),keys.end(),keyword[ESURF].c_str());
  _esurf = atof((*iter).value());
  // nfreeze 
  iter = find(keys.begin(),keys.end(),keyword[NFREEZE].c_str());
  _nfreeze = atoi((*iter).value().c_str()); 
  // iwrite 
  iter = find(keys.begin(),keys.end(),keyword[WRITE].c_str());
  _iwrite = atoi((*iter).value()); 
  // coordname 
  if(iter != keys.end()) keys.erase(iter);  
  // restartname 
  iter = find(keys.begin(),keys.end(),keyword[RESTART]);
  if(iter != keys.end()) keys.erase(iter);
  // trajecotry
  iter = find(keys.begin(),keys.end(),keyword[TRAJ]);
  if(iter != keys.end()) keys.erase(iter);
  // energy
  iter = find(keys.begin(),keys.end(),keyword[ENERGY]);
  if(iter != keys.end()) keys.erase(iter);
  // log_file
  iter = find(keys.begin(),keys.end(),keyword[LOG]);
  if(iter != keys.end()) keys.erase(iter);
  // err_file
  iter = find(keys.begin(),keys.end(),keyword[ERR]);
  if(iter != keys.end()) keys.erase(iter);
  // initial
  iter = find(keys.begin(),keys.end(),keyword[CALC]);
  if((*iter).value()== keyword[INITRUN].c_str()) _calctype = 0;
  else if((*iter).value() == keyword[CONTRUN].c_str()) _calctype = 1;
  else {
    sprintf(line,"Error in input file \"%s\" calculation type = \"%s\"",
	    _inname.c_str(),(*iter).value()); writeErr(line);
    sprintf(line,"Can be either \"%s\" or \"%s\"\n",
	    keyword[INITRUN].c_str(),keyword[CONTRUN].c_str()); 
    writeErr(line);
    exit(1);
  }
  // zerotm
  iter = find(keys.begin(),keys.end(),keyword[TMZERO]);
  if(iter != keys.end()) keys.erase(iter);
  
  // cutoff
  iter = find(keys.begin(),keys.end(),keyword[CUTOFF]);
  if(iter != keys.end()) keys.erase(iter);

  // paramfile
  iter = find(keys.begin(),keys.end(),keyword[PARAM]);
  if(iter != keys.end()) keys.erase(iter);

  // systemfile
  iter = find(keys.begin(),keys.end(),keyword[SYSTEM]);
  if(iter != keys.end()) keys.erase(iter);
  
  // for(iter = keys.begin();iter!=keys.end();++iter) cout << *iter <<endl;
  
#endif  

  setStderr(errname().c_str());
  this->print(setLogfile(logname().c_str()));

  // a few consistancy checks
  if(_rcut > _xbox/2){
    sprintf(line,"WARNING: Rut off distance (%g) > (%g) 1/2 Box size!!",
	    _rcut, _xbox/2);
    writeErr(line);
  } 
  if(_nfreeze >= _nmol){ 
    sprintf(line,"ERROR: number of frozen atoms %d dosn't compute! %d>=%d!",
	    _nfreeze,_nfreeze,_nmol);
    writeErr(line);
    exit(1);
  }
  _dof = DIM*(_nmol);
}
