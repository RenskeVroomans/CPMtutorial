#ifndef MiscHeader
#define MiscHeader


#include "Header.hh"

double returnNormal(double mean, double stdev);

double Derivative(double ht,double enhancer,double decay, double x);
double RungeKutta2(double HT, double enhancer, double decay, double x);
double RungeKutta4(double HT, double enhancer, double decay, double x);

#endif
