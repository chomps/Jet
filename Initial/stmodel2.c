
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

double f( double r , double a , double b , double M ){
   return( (M/M_PI/a)/( (r-b)*(r-b) + 1./a/a ) );
}

void initial( double * prim , double * x ){

   double Msun = 1.0;
   double R0   = 1.0;
   double r  = x[0];

   double rhoc = 3e7*Msun/pow(R0,3.);

   double R1 = .0017*R0;
   double R2 = .0125*R0;
   double R3 = .65*R0;

   double x1 = 0.00354*R0;
   double x2 = 0.01385*R0;
   double x3 = 0.112*R0;

   double k1 = 3.24;
   double k2 = 2.57;
   double n  = 16.7;

   double rho_wind = 1e-9*Msun/pow(R0,3.);

   double X = 1.-r/R3;
   if( X<0.0 ) X = 0.0;

   double rho = rhoc*pow(X,n)/( 1. + pow(r/R1,k1)/(1.+pow(r/R2,k2)) ) + rho_wind*pow(R3/r,2.);

   double Pp = 1e-6*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = 0.0;

}
