
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

double f( double r , double a , double b , double M ){
   return( (M/M_PI/a)/( (r-b)*(r-b) + 1./a/a ) );
}

void initial( double * prim , double * x ){

   double w = 1.3;

   double v_max = 0.1;

   double M0 = 1.0;
   double R0 = 1.0;
   double r  = x[0];
   double th = x[1];
   
   r *= pow(w,2./3.)*sqrt( cos(th)*cos(th) + sin(th)*sin(th)/w/w );

   double R1 = .06*R0;
   double rhoc = M0/pow(R0,3.)/.0102913;  ///.003479;

//	2	.0221136
//	2.5	.0102913
//	3	.00554687
//	3.5	.00347901

   double rhostar = rhoc/( 1. + pow(r/R1,2.5) )*pow(1.-r/R0,.75);
   if( r>R0 ) rhostar = 0.0;;

   double rho_trans = 1e-8;
   double rho_ISM = 1e-20;
   double r_ext = x[0];

   double rho_ext = rho_trans*exp(-r_ext/R0) + rho_ISM;
   double rho = rhostar + rho_ext;

   double v = v_max*(x[0]/R0)*rhostar/rho;
   double u = v/sqrt(1.-v*v);

   double Pp = 1e-6*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = u;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = 0.0;

}
