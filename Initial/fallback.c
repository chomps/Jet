
#include "../paul.h"

static double k = 1.0;

void setICparams( struct domain * theDomain ){
   k = theDomain->theParList.Explosion_Energy;
}

void initial( double * prim , double * x ){

   double r  = x[0];

   double M = 1.0;
   double E = 1.0;
   double R = 1.0;

   double klog = k*log(r/R);

   double rho = (k/M_PI)*M/(4.*M_PI*r*r*r)/( 1. + klog*klog ); 
   double P0 = 1e-3*rho;
   double Rexp = 1e-4*R;
   double Pexp = E/(Rexp*Rexp*Rexp)*pow( 2.*M_PI , -1.5 )/1.5;
   double Pp = P0 + Pexp*exp(-0.5*r*r/(Rexp*Rexp) );

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ){
      double X = 0.0;
      if( r>R ) X = 1.0;
      prim[NUM_C] = X;
   }
}
