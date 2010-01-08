/*! \file simkeys.cc
  file holds the public functions for the simkey class
*/

#include "simkeys.h++"

//using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

std::ostream & operator<<(std::ostream& stream,SIMKEY & src)
{
  stream << "key:"<<src.isReq()<<", "<<src.ktype()<<", "<<src.iset()<<", \""
	 <<src.key() <<"\", \""<<src.value()<<"\", \""<<src.defvalue()<<"\"";
  return stream;
}

void SIMKEY::print(FILE *fp=stdout)
{
  fprintf(fp,"key: %d %d %d \"%s\", \"%s\", \"%s\"\n",
	  _isReq,_ktype,_iset,_key.c_str(),_value.c_str(),_defvalue.c_str());
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

SIMKEY & SIMKEY :: operator=(const SIMKEY &src)
{
  // so we don't assign to ourselves and use subclasses not to memory leak
  if(this == &src) return (*this);
  
  _isReq=src._isReq;
  _ktype=src._ktype;
  _key = src._key;
  _value = src._value;
  _defvalue = src._defvalue;
  return (*this);
}

SIMKEY :: SIMKEY(void) //!< constructor
{ 
  _iset=_isReq=0;
  _ktype=SIMKEY::NONE;
  // _key=_value=_defvalue=NULL;
}

SIMKEY :: ~SIMKEY(void) //!< destructor
{
}

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

  _ktype=tp;
  _key = ky;
  _value = dv;
  _defvalue = dv;
}
