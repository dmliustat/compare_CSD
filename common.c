#define _USE_MATH_DEFINES
#include "header.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


static double const normalization = 0.5 * M_SQRT1_2 * M_2_SQRTPI;

double *dou_vec_init(const size_t num);
long *long_vec_init(const size_t num);
short *sh_vec_init(const size_t num);

unsigned long seed[2];


extern long L;

double *dou_vec_init(const size_t num)
{
  double *p;
  if((p = (double *) calloc(num, sizeof(double))) == NULL)
    perror("initialise_array():calloc failed!");
  return(p);
}


long *long_vec_init(const size_t num)
{
  long *p;
  if((p = (long *) calloc(num, sizeof(long))) == NULL)
    perror("initialise_array():calloc failed!");
  return(p);
}


short *sh_vec_init(const size_t num)
{
  short *p;
  if((p = (short *) calloc(num, sizeof(short))) == NULL)
    perror("initialise_array():calloc failed!");
  return(p);
}


