
#include "../paul.h"

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];
   double z = r*cos(th);

   double H  = .189;
   double a  = .299;
   double b  = .001;
   double z0 = 1.0;
   double kz = 1.0;

   double rho = a*exp(-.5*z*z/H/H) + b*pow( 2.*z0/(fabs(z)+z0), kz );

   double rin  = R_MIN;
   double Pmax = 0.00000625/pow(rin,3.);
   double Pp;

   if( r<rin ) Pp = Pmax;
   else Pp = 1e-8*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}
