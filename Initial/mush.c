
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

   //The following is consistent with an initial
   //time t = r0/vmax = .05477

   double r = x[0];
   double r0 = 0.5;
   double rho0 = 1.0;
   double Pmin = 1e-5;

   double rho = rho0;
   double X = 0.0;
   double v = 0.0;

   if( r < r0 ){
      rho = .5;
      X = 1.0;
   }
   double r_max = 1.0;
   double P0 = 1.0;
   double Pp = Pmin + (P0-Pmin)*.5*(cos(M_PI*r/r_max)+1.);
 
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = v;
   prim[UU2] = 0.0;
   if( NUM_N > 0 ) prim[NUM_C] = X;

}
