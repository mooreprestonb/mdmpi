
/*! \file mdinclude.h++
  \brief file that contains useful macros and constants and defines

 */
#ifndef _MDINCLUDE_H
#define _MDINCLUDE_H


/*! \def anint(A,B)
    \brief A macro that returns the nearest integer of \a A.
*/
#ifndef anint
#define anint(A) ((int)((A)+((A)>=0?.5:-.5)))
#endif

/*! \def MAX(A,B)
    \brief A macro that returns the maximum of \a A and \a B.
*/
#define MAX(A,B) (((A)>(B))?(A):(B))

/*! \def MIN(A,B)
    \brief A macro that returns the mimimum of \a A and \a B.
*/
#define MIN(A,B) (((A)<(B))?(A):(B))

const int DIM = 3;
const int WDIM = 3;

#define MAXALLOC 1000

#define NCONF 10
// #define NCONF 1
// #define DEBUG
#define WRITE_BOX
#define ZEROTM
#define RESETPOS

const int MAXLINE=256;

const double DX=1.e-5;

#endif
