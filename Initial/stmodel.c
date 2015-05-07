
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

double f( double r , double a , double b , double M ){
   return( (M/M_PI/a)/( (r-b)*(r-b) + 1./a/a ) );
}

void initial( double * prim , double * x ){

   double M  = 1.0;
   double R  = 1.0;
   double r  = x[0];

   double rho0 = 1e-6;

   double a1 = 23.4*(0.4/R);
   double b1 = 0.058*(R/0.4);
   double M1 = 11.27/11.2*M;

   double a2 = 472.*(0.4/R);
   double b2 = 0.002*(R/0.4);
   double M2 = 3.24/11.2*M;

   double f1 = f( r , a1 , b1 , M1 ); 
   double f2 = f( r , a2 , b2 , M2 ); 

   //double Vol = 2.*pow( sqrt(2.*M_PI)*R , 3. );
   //double rho0 = M/Vol;

   //double rho = rho0*exp(-.5*r*r/R/R);

   double rho = (f1+f2)/4./M_PI/r/r * pow( 1. - r/R , 3.85 ) + rho0;

   if( r >= R ) rho = rho0*(R*R/r/r);

   double Pp = 1e-4*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = 0.0;

}
