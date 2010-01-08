/*! \file simkeys.h++
  file defines the key type...
*/

#ifndef _SIMKEYS_H
#define _SIMKEYS_H

#include <string>
#include <iostream>
#include <stdio.h>

class SIMKEY
{
private:
  int _isReq,_ktype,_iset;
  std::string _key,_value,_defvalue;
public:
  enum {NONE,INT,DOUBLE,STRING};
  SIMKEY(void); // constructor
  SIMKEY (const SIMKEY &src){*this = src;} // copy constructor
  ~SIMKEY(void); // destructor
  void setdef(int,int,const std::string &,const std::string &);
  int isReq(void){ return _isReq;}
  int iset(void){ return _iset;}
  int ktype(void){ return _ktype;}
  void key(const std::string &key){ _key = key;}
  void set(const std::string &,const std::string &);
  std::string key(void){return _key;}
  void value(const std::string &value){ _value = value;}
  std::string value(void){return _value;}
  void defvalue(const std::string &defvalue){ _defvalue = defvalue;}
  std::string defvalue(void){return _defvalue;}
  SIMKEY & operator=(const SIMKEY &src); // copy assignment
  SIMKEY & operator+=(const SIMKEY &src); // add value
  bool operator==(const SIMKEY &);
  bool operator!=(const SIMKEY &);
  bool operator!=(const std::string &);
  bool operator==(const std::string &);
  void print(FILE *);
};

std::ostream& operator<<(std::ostream &,SIMKEY &);// output

#endif
