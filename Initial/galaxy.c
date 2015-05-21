
#include "../paul.h"

static double R_MIN = 0.0;

void setICparams( struct domain * theDomain ){
   R_MIN = theDomain->theParList.rmin;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];
   double z = r*cos(th);

   double H  = .189;
   double a  = .299;
   double b  = .001;
   double z0 = 1.0;
   double kz = 2.0;

   double rho = a*exp(-.5*z*z/H/H) + b*pow( z0/(fabs(z)+z0), kz );

   double rin  = 1.*R_MIN;
   double Pmax = 1.0/pow(rin,3.);
   double Pp = 1e-8*rho;

   Pp += Pmax*exp(-.5*(r*r/rin/rin));

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}
