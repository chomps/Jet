
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

double f( double r , double a , double b , double M ){
   return( (M/M_PI/a)/( (r-b)*(r-b) + 1./a/a ) );
}

void initial( double * prim , double * x ){

   double M0 = 1.0;
   double R0 = 1.0;
   double r  = x[0];

   //double rhoc = 0.076415*M0/pow(R0,3.);
   double rhoc = 0.085*M0/pow(R0,3.);

   double k = 6.0;
   double k2 = 1.5;
   double k3 = 2.5;

   double n = 0.5;
   double A = 0.00583*M0/pow(R0,3.);

   double Rmax = 1.0*R0;
   
   double rho_wind = 1e-9*M0/pow(R0,3.);

   double X = 1.-r/Rmax;
   if( X<0.0 ) X = 0.0;

   double rho = rhoc*pow(r/R0,-k3)*pow(X,k2)/pow( 1. + pow(A/rhoc*pow(r/R0,-k) ,-n) ,1./n) + rho_wind*pow(Rmax/r,2.);

   double rho_max = 1e6;
   if( rho > rho_max ) rho = rho_max;

   double Pp = 1e-4*rho;
   if( Pp > 10.*M0/R0/R0/R0 ) Pp = 10.*M0/R0/R0/R0;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = 0.0;
   if( NUM_N > 1 ){
      if( r > 0.03 ) prim[NUM_C+1] = 0.0;
      else prim[NUM_C+1] = 1.0;
   }

}
