/*! \file simkeys.h++
  file defines the key type...
*/

#ifndef _SIMKEYS_H
#define _SIMKEYS_H

// #define NOSTRING

#ifdef NOSTRING
#else
#include <string>
#endif

#include <iostream>
#include <stdio.h>

class SIMKEY
{
private:
  int _isReq,_ktype,_iset;
#ifdef NOSTRING
  char *_key,*_value,*_defvalue;
#else
  std::string _key,_value,_defvalue;
#endif
  void freeValue(void);
public:
  enum {NONE,INT,DOUBLE,STRING};
  SIMKEY(void); // constructor
  SIMKEY (const SIMKEY &src){*this = src;} // copy constructor
  ~SIMKEY(void); // destructor
  int isReq(void){ return _isReq;}
  void isReq(const int ireq){ _isReq = ireq;}
  int iset(void){ return _iset;}
  void iset(const int iset){_iset = iset;}
  int ktype(void){ return _ktype;}
#ifdef NOSTRING
  void setdef(int,int,const char *,const char *);
  void key(const char *key){_key = strdup(key);}
  void set(const char *,const char *);
  char *key(void){return _key;}
  void value(const char *value)
  {
    if(_value != NULL) free(_value);
    _value = strdup(value);
  }
  char *value(void){return _value;}

  void defvalue(const char *defvalue)
  {
    if(_defvalue != NULL) free(_defvalue);
    _defvalue = strdup(defvalue);
  }
  char *defvalue(void){return _defvalue;}
  
  bool operator!=(char const *);
  bool operator==(char const *);
#else
  void setdef(int,int,const std::string &,const std::string &);
  void key(const std::string &key){ _key = key;}
  void set(const std::string &,const std::string &);
  std::string key(void){return _key;}
  void value(const std::string &value){ _value = value;}
  std::string value(void){return _value;}
  void defvalue(const std::string &defvalue){ _defvalue = defvalue;}
  std::string defvalue(void){return _defvalue;}
  bool operator!=(const std::string &);
  bool operator==(const std::string &);
#endif

  SIMKEY & operator=(const SIMKEY &src); // copy assignment
  SIMKEY & operator+=(const SIMKEY &src); // add value
  bool operator==(const SIMKEY &);
  bool operator!=(const SIMKEY &);
  void print(FILE *);
};

std::ostream& operator<<(std::ostream &,SIMKEY &);// output

#endif
