
#include "../paul.h"

static double t = 0.0;

void setICparams( struct domain * theDomain ){
   t = theDomain->theParList.t_min;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double E_ej = 0.01;
   double M_ej = 1.0;

   double n = 10.0;
   double d = 2.0;

//Ejecta Density
   double vt2 = E_ej/M_ej*2.*(5.-d)*(n-5.)/((3.-d)*(n-3.));
   double vt = sqrt(vt2);
   double R = vt*t;

   double zeta = ( (n-3.)*(3.-d)/(n-d) )/4./M_PI;
   double rho_T = zeta*M_ej/pow(R,3.);

   double rho_ej = rho_T*pow(r/R,-d);
   if( r>R ) rho_ej = rho_T*pow(r/R,-n);

   double rho = rho_ej;
   double Pp = 1e-5*rho*vt2;
   double v = r/t;
   if( v > 0.5 ) v = 0.0;
   double u = v/sqrt(1.-v*v);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = u;
   prim[UU2] = 0.0;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = 0.0;
   }
   if( NUM_N>0 && r<R ) prim[NUM_C] = 1.0;

}


