/*! \file file to parse flags (ie command line arguments) 

 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "flags.H"

/*! help the user out by showing the command line arguements */
static void usage(char s[])
{
  fprintf(stderr,"USAGE: %s [-hq] <infile>\n",s);
  exit(1);
}

/*! actually parse up the command line */
void FLAGS :: read(int argc,char **argv)
{
  /* set default flags */
  _argv = argv;
  _argc = argc;
  _command = strdup(argv[0]);

  while(argc > 1){
    if((*++argv)[0] == '-'){
      --argc;
      while(*(++argv[0]) != (char)NULL ){
        switch(*(argv[0])){
	case 'q':
	  _flq = 1;
	  break;
//         case 'n':
// 	  sscanf(*(++argv),"%d",&_npt);
// 	  --argc; *(++argv[0]) = (char)NULL;--argv[0];
//           break;
//         case 'x':
// 	  sscanf(*(++argv),"%lg",&_xmax);
// 	  --argc; *(++argv[0]) = (char)NULL;--argv[0];
//           break;
        case 'h':
          fprintf(stderr,"Program %s run molecular dynanics code",_command);
          fprintf(stderr,"from input file\n");
          fprintf(stderr,"-q = quite\n");
          fprintf(stderr,"-h = help (this screen)\n");
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
    } else {
      usage(_command);
    }
  }
  if(_infile==NULL){
    usage(_command);
  }
}

