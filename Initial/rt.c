
#include "../paul.h"

static double R_MIN = 0.0;
static double GAMMA_LAW = 0.0;
static double ETA = 0.0;

void setICparams( struct domain * theDomain ){
   R_MIN = theDomain->theParList.rmin;
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
   ETA = theDomain->theParList.Explosion_Energy;
   srand(666);
   rand();
}

void initial( double * prim , double * x ){

   //The following is consistent with an initial
   //time t = r0/vmax = .05477

   double E = 1.0;
   double M = 1.0;

   double r = x[0];
   double r0 = 0.01;
   double rho0 = 1.0;
   double Pmin = 1e-5;

   double rho = rho0;
   double v   = 0.0;
   double X   = 0.0;

   double V = 4./3.*M_PI*r0*r0*r0;
   double vmax = sqrt(10./3.*E/M);

   if( r < r0 ){
      rho = M/V;
      v = vmax*r/r0;
      X = 1.0;
   }
   double vpert = 1e-3*v*((double)rand()/(double)RAND_MAX-.5);
 
   prim[RHO] = rho;
   prim[PPP] = Pmin;
   prim[UU1] = v + vpert;
   prim[UU2] = 0.0;
   if( NUM_C > 4 ) prim[UU3] = 0.0;
   if( NUM_N > 0 ) prim[NUM_C] = X;

}
