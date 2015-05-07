
#include "../paul.h"

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double eta = 10.0;
   double w0 = 0.1;
   double n = 20.0;

   double rm = .0003;
   double r0 = .1;

   double rho0 = 1.0;
   double M = rho0*(4.*M_PI/3.)*r0*r0*r0*r0/rm;
   double E = eta*M;

   double rho = rho0*( pow( r0/(r+rm) , 4. ) + 1. );
   double e0 = E/pow(sqrt(M_PI)*rm,3.);
   double P0 = (GAMMA_LAW-1.)*e0;
   double Pmin = 1e-5*rho;

   double Pp = P0*exp(-r*r/rm/rm) + Pmin;

   double X = 4.*(r-r0)/r0;
   double vr = w0*exp(-pow(X,4.))*sin(n*th);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = vr;
   prim[UU2] = 0.0;

}
