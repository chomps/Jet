
#include "../paul.h"

/*
Analytic parametric solution for the density and radius:
(t goes from small to 1.0)

time = t*sqrt(1./t-1.) 
+ .5*atan(sqrt(1./t-1.)*(2.*t-1.)/2./(t-1.))
+ pi/4.

all over sqrt(8*pi*G/3.)

rho   = t**-3
r_out = t
*/

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r = x[0];
   double R = 1.0;
   double rho1 = 1.0;
   double rho0 = 1e-3;
   double P0 = 1e-9;

   double rho = rho0;
   if( r<R ) rho = rho1;

   prim[RHO] = rho;
   prim[PPP] = P0;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}
