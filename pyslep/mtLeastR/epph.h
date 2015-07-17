#ifndef _EPPH_H_
#define _EPPH_H_

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define delta 1e-8

#define innerIter 1000
#define outerIter 1000


/*
  This is a head file for the used C files

*/





void eplb(double * x, double *root, int * steps, double * v,int n, double z, double lambda0);

void  epp1(double *x, double *v, int n, double rho);

void  epp2(double *x, double *v, int n, double rho);

void  eppInf(double *x, double * c, int * iter_step, double *v,  int n, double rho, double c0);

void zerofind(double *root, int * iterStep, double v, double p, double c, double x0);

double norm(double * v, double p, int n);

void  eppO(double *x, double * cc, int * iter_step, double *v,  int n, double rho, double p);

void epp(double *x, double * c, int * iter_step, double * v, int n, double rho, double p, double c0);

#endif // _EPPH_H_
