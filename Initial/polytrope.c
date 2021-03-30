
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r = x[0];

   double R = 1.0;
   double G = 1.0;
   double M = 1.0;
   double rhoc = 1.0;

   double rho = 1.0;
   double Pp  = 1.0;

   double rho_min = 1e-8;

   int n=5;

   if( n==1 ){
      rhoc = M/R/R/R/4./M_PI/M_PI;
      rho = rhoc*sin(r/R)/(r/R) + rho_min;
      if( r > M_PI*R ) rho = rho_min;
      Pp = 2.*M_PI*G*rho*rho*R*R;
   }else if (n==5){
      rhoc = M/R/R/R/4./M_PI/sqrt(3.);
      rho = rhoc/pow( 1. + r*r/R/R/3. , 2.5 );
      Pp = 2./3.*M_PI*G*rhoc*rhoc*R*R/pow( 1. + r*r/R/R/3. , 3. );
   }

   prim[RHO] = rho;
   prim[PPP] = 1.065*Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}


