
#include "../paul.h"

static double R_MIN = 0.0;
static double GAMMA_LAW = 0.0;
static double ETA = 0.0;

void setICparams( struct domain * theDomain ){
   R_MIN = theDomain->theParList.rmin;
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
   ETA = theDomain->theParList.Explosion_Energy;
}

void initial( double * prim , double * x ){

   double eta = ETA;

   double r = x[0];
   double r0 = 100.*R_MIN;
   double rho0 = 1.0;
   double Pmin = 1e-5;

   double E = 1.0;
   double e0 = E/pow(sqrt(M_PI)*r0,3.);

   double ex   = e0*exp(-r*r/r0/r0);
   double Px   = (GAMMA_LAW-1.)*ex;
   double rhox = ex/eta;

   double rho,Pp;
   rho = rhox + rho0;
   Pp = Px + Pmin;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}
