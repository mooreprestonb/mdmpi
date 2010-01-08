/*! \file simkeys.cc
  file holds the public functions for the simkey class
*/

#include "simkeys.h++"

//using namespace std;

#ifdef NOSTRING
#else
#include <string>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

void SIMKEY :: freeValue(void)
{
#ifdef NOSTRING
  if(_key!=NULL) free(_key);
  if(_value!=NULL) free(_value);
  if(_defvalue!=NULL) free(_defvalue);
  _key=_value=_defvalue=NULL;
#endif
}

std::ostream & operator<<(std::ostream& stream,SIMKEY & src)
{
  stream << "key: req="<<src.isReq()<<", type="<<src.ktype()<<", set="<<src.iset()<<", \""
#ifdef NOSTRING
	 <<src.key() <<"\", \""<<src.value()<<"\", \""<<src.defvalue()<<"\"";
#else
	 <<src.key() <<"\", \""<<src.value()<<"\", \""<<src.defvalue()<<"\"";
#endif
  return stream;
}

void SIMKEY::print(FILE *fp=stdout)
{
  fprintf(fp,"key: req:%d type:%d set:%d ",_isReq,_ktype,_iset);
#ifdef NOSTRING
  fprintf(fp,"\"%s\", \"%s\", \"%s\"\n",_key,_value,_defvalue);
#else
  fprintf(fp,"\"%s\", \"%s\", \"%s\"\n",_key.c_str(),_value.c_str(),_defvalue.c_str());
#endif
}

bool SIMKEY::operator==(const SIMKEY &src)
{
  if(this == &src) return true; // same address 
  if(_isReq != src._isReq) return false;
  if(_ktype != src._ktype ) return false;
  if(_key != src._key) return false;
  if(_value != src._value) return false;
  if(_defvalue!= src._defvalue) return false;

  return true;
}

SIMKEY & SIMKEY::operator+=(const SIMKEY &src1)
{
  if(this != &src1){ 
    if(_key == src1._key){
      _iset++; 
      _value = src1._value;
    } else fprintf(stderr,"ERROR: adding different keys?\n");
  } else { // adding ourselves ie a += a;
    fprintf(stderr,"ERROR: adding to ourselves in += op of simkey?\n");
  }
  return (*this);
}

bool SIMKEY::operator!=(const SIMKEY &src1)
{
  if(*this == src1) return false;
  return true;
}

#ifdef NOSTRING
bool SIMKEY::operator==(char * const tkey)
{
  if(strcasecmp(tkey,_key)==0) return true;
  return false;
}
bool SIMKEY::operator!=(char * const tkey)
{
  return !(*this==tkey);
}
#else
bool SIMKEY::operator==(const std::string &tkey)
{
  if(tkey == _key) return true;
  return false;
}

bool SIMKEY::operator!=(const std::string &tkey)
{
  if(tkey == _key) return false;
  return true;
}
#endif

SIMKEY & SIMKEY :: operator=(const SIMKEY &src)
{
  // so we don't assign to ourselves and use subclasses not to memory leak
  if(this == &src) return (*this);
  
  _isReq=src._isReq;
  _ktype=src._ktype;
  _iset =src._iset;
#ifdef NOSTRING
  this->freeValue();
  _key = strdup(src._key);
  _value = strdup(src._value);
  _defvalue = strdup(src._defvalue);
#else
  _key = src._key;
  _value = src._value;
  _defvalue = src._defvalue;
#endif
  return (*this);
}

SIMKEY :: SIMKEY(void) //!< constructor
{ 
  _iset=_isReq=0;
  _ktype=SIMKEY::NONE;
#ifdef NOSTRING
  _key=_value=_defvalue=NULL;
#endif
}

SIMKEY :: ~SIMKEY(void) //!< destructor
{
#ifdef NOSTRING
  this->freeValue();
#endif
}

#ifdef NOSTRING
//! sets the key and value
void SIMKEY::set(char *ky,char *va)
{
  if(_key!=NULL) free(_key);
  if(_value!=NULL) free(_value);
  _key = strdup(ky);
  _value = strdup(va);
}

void SIMKEY::setdef(int req,int tp,const char *ky,const char *dv)
{
  if(req!=0 && req!=1){
    fprintf(stderr,"ERROR in SIMKEY::setdef var::req is not 0 or 1 but %d\n",
	    req);
    req=0;
  }
  _isReq=req;
  if(tp!=NONE && tp!=INT && tp!=DOUBLE &&tp!=STRING){
    fprintf(stderr,"ERROR in SIMKEY::setdef type variable is a valid type\n");
    fprintf(stderr,"But is %d\n",tp);
    tp=NONE;
  }

  _iset=0;
  _ktype=tp;

  this->freeValue();
  if(ky != NULL) _key = strdup(ky);
  if(dv != NULL) _value = strdup(dv);
  if(dv != NULL) _defvalue = strdup(dv);
}

#else
//! sets the key and value
void SIMKEY::set(const std::string &ky,const std::string &va)
{
  _key = ky;
  _value = va;
}

void SIMKEY::setdef(int req,int tp,
		    const std::string & ky,const std::string & dv)
{
  if(req!=0 && req!=1){
    fprintf(stderr,"ERROR in SIMKEY::setdef var::req is not 0 or 1 but %d\n",
	    req);
    req=0;
  }
  _isReq=req;
  if(tp!=NONE && tp!=INT && tp!=DOUBLE &&tp!=STRING){
    fprintf(stderr,"ERROR in SIMKEY::setdef type variable is a valid type\n");
    fprintf(stderr,"But is %d\n",tp);
    tp=NONE;
  }

  _iset=0;
  _ktype=tp;
  _key = ky;
  _value = dv;
  _defvalue = dv;
}

#endif
