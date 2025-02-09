#include "Misc.hh"

double returnNormal(double mean, double stdev)
{
  static int spareavailable=0;
  static double spare=0.;
  double val1, val2;
  double s, fac;
  if(spareavailable){
    spareavailable=0;
    return mean +spare*stdev;
  }
  else{
    do{
      
      val1=uniform()*2. -1.;
      val2=uniform()*2. -1.;
      s=val1*val1+val2*val2;
      
    }while (s>=1.||s==0);
    
    fac=sqrt((-2.*log(s))/s);
    
    spare=val1*fac;
    spareavailable=1;
    return mean+val2*fac*stdev;

  }
} 

double Derivative(double ht,double enhancer,double decay, double x)
{
  double deriv;
  double y;
 
  deriv=enhancer-decay*x;

  y=ht*deriv;
  return y;
}

double RungeKutta2(double HT,double enhancer, double decay, double x)
{
  double k1;
  double k2;
  double rk;

  k1=Derivative(HT,enhancer,decay,x);
  k2=Derivative(HT,enhancer,decay,x+k1);

  rk=0.5*(k1+k2);
  return rk;
}


double RungeKutta4(double HT, double enhancer, double decay, double x)
{
  double k1;
  double k2;
  double k3;
  double k4;
  double rk;

  k1=Derivative(HT,enhancer,decay,x);
  k2=Derivative(HT,enhancer,decay,x+0.5*k1);
  k3=Derivative(HT,enhancer,decay,x+0.5*k2);
  k4=Derivative(HT,enhancer,decay,x+k3);
  rk=(1./6.)*(k1+2*k2+2*k3+k4);

  return rk;
}