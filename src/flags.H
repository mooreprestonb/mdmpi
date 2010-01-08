/*! \file flags.H
  This file holds a simple class to parse the command lines arguments
  and io functions to return command line flags
*/

#ifndef _FLAGS_H
#define _FLAGS_H

//! \class FLAGS store and sort command line arguements 
class FLAGS {
private:
  //! command line arguments
  char ** _argv;
  //! # of command line arguments
  int _argc;
  //! command name
  char * _command;
  //! input file name
  char * _infile;
  //! quite/verbose flag
  int _flq;
public:
  //! construct by zero everything
  FLAGS(void) {_command=_infile=NULL;_flq= 0; }
  //! destruct which frees strings if necessary
  ~FLAGS(void) {
    if(_command != NULL) free(_command);
    if(_infile != NULL) free(_infile);
  }
  //! read the command lines
  void read(int, char **);
  //! return pointer to the input file
  char* infile(void){return _infile;}
  //! return pointer to the command name
  char* command(void){return _command;}
  //! return the quite flag
  int q(void){return _flq;}
  //! return the number of command line arguements
  int argc(void){return _argc;}
  //! return command line arguements
  char ** argv(void){return _argv;}
};

#endif
