#ifndef GeneralHeader
#define GeneralHeader


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dSFMT.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <cmath>
#include <list> 
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <algorithm>
#include <ext/numeric>
#include <boost/utility.hpp>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdlib.h> 
#include <cctype>
#include <sstream>

#define TRUE 1
#define FALSE 0


#define RUN
//#define INITGENOME

#define SUMINTEGRATION

typedef struct POSITION{
  int xx;
  int yy;
}Position;

typedef struct FPOSITION{
  double xx;
  double yy;
}FPosition;

using namespace std;
using namespace __gnu_cxx;


extern char despath[500];

extern int seed;

extern dsfmt_t dsfmt;
inline double uniform() { return dsfmt_genrand_close_open(&dsfmt); }

const int randlength=10000000;
extern double randarray[randlength];
extern int randloc;
extern int startloc[10];


#endif


