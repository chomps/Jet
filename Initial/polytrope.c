
#include "../paul.h"

static double f_le[100];

void setICparams( struct domain * theDomain ){

   int i;
   for( i=0 ; i<100 ; ++i ){
       //Solve Lane Emden for f(x)
   }

}

void initial( double * prim , double * x ){

   double x = 6.3*r/R;
   int nx = (int)( x*100 );

   double rho = rho0*f_le[nx];

   prim[RHO] = rho;
   prim[PPP] = rho*1e-6;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}


