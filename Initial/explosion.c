
#include "../paul.h"

static double R_MIN = 0.0;
static double GAMMA_LAW = 0.0;

void setICparams( struct domain * theDomain ){
   R_MIN = theDomain->theParList.rmin;
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double r0 = 0.02;
   double rho = 1.0;
   double Pmin = 1e-6;

   double E = 1.0;
   double Px   = (GAMMA_LAW-1.)*E*exp(-r*r/r0/r0)/pow(sqrt(M_PI)*r0,3.);

   double Pp = Pmin + Px;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}
