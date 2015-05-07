
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   double r = x[0];
   double rho,Pp;
   if( r < 0.25 ){
      rho = 1.0;
      Pp  = 1.0;
   }else{
      rho = 0.1;
      Pp  = 0.1;
   }
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
}
