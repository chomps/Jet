
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r = x[0];

   double r0 = 0.1; 
   double T0 = 2.5; 

   double rho_ej = pow(r0,-3.)*exp(-pow(r/r0,2.));
   double rho_ism = 1.0; 
   double rho = rho_ej + rho_ism;

   double X = rho_ej/rho;   
   double v = r/r0*X;


   prim[RHO] = rho;
   prim[PPP] = rho*T0;
   prim[UU1] = v;
   prim[UU2] = 0.0;
   if( NUM_N > 0 ) prim[NUM_C] = X;

}
