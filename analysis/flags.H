
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef _FLAGS_H
#define _FLAGS_H

const int NPT=100;
const double XMIN=0.;
const double XMAX=10.;

class FLAGS {
private:
  char *_command,*_infile,*_outfile;
  int _flq,_npt;
  double _xmin,_xmax;
public:
  FLAGS(void) {
    _command=_infile=_outfile=NULL;
    _flq= 0;_npt=NPT;_xmin = XMIN;_xmax = XMAX;
  }
  ~FLAGS(void) {
    if(_command != NULL) free(_command);
    if(_infile != NULL) free(_infile);
    if(_outfile != NULL) free(_outfile);
  }
  void read(int, char **);
  char* infile(void){return _infile;}
  char* outfile(void){return _outfile;}
  char* command(void){return _command;}
  int n(void){return _npt;}
  int q(void){return _flq;}
  double xmin(void){return _xmin;}
  double xmax(void){return _xmax;}
};

/*-----------------------------------------------------------------------*/
void usage(char s[])
{
  fprintf(stderr,"USAGE: %s [-hq] [-n #] [-x min] <infile> <outfile>\n",s);
  exit(1);
}

/*-----------------------------------------------------------------------*/
void FLAGS :: read(int argc,char *argv[])
{
  /* set default flags */
  _command = strdup(argv[0]);

  while(argc > 1){
    if((*++argv)[0] == '-'){
      --argc;
      while(*(++argv[0]) != (char)NULL ){
        switch(*(argv[0])){
	case 'q':
	  _flq = 1;
	  break;
        case 'n':
	  sscanf(*(++argv),"%d",&_npt);
	  --argc; *(++argv[0]) = (char)NULL;--argv[0];
          break;
        case 'x':
	  sscanf(*(++argv),"%lg",&_xmax);
	  --argc; *(++argv[0]) = (char)NULL;--argv[0];
          break;
        case 'm':
	  sscanf(*(++argv),"%lg",&_xmin);
	  --argc; *(++argv[0]) = (char)NULL;--argv[0];
          break;
        case 'h':
          fprintf(stderr,"Program %s generates a g(r) ",_command);
	  fprintf(stderr,"of x y z data which reads from stdin or infile\n");
          fprintf(stderr,"outputs to stdout or outfile\n");
          fprintf(stderr,"-q = quite\n");
          fprintf(stderr,"-h = help (this screen)\n");
          fprintf(stderr,"-n # = number of points to use\n");
          fprintf(stderr,"-x # = max of x\n");
          fprintf(stderr,"-m # = min of x\n");
          usage(_command);
          break;
        default:
          usage(_command);
          break;
        }
      }
    } else if (_infile == NULL){
      _infile = strdup(argv[0]);
      --argc;
    } else if (_outfile == NULL) {
      _outfile = strdup(argv[0]);
      --argc;
    } else {
      usage(_command);
    }
  }
  if(!_flq){
    fprintf(stdout,"infile = \"%s\"\n",_infile);
    fprintf(stdout,"outfiles = \"%s\"\n",_outfile);
    fprintf(stdout,"npt = %d\n",_npt);
    fprintf(stdout,"max = %lg\n",_xmax);
    fprintf(stdout,"min = %lg\n",_xmin);
  }
}

#endif
