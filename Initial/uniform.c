
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   prim[RHO] = 1.0e-7;
   prim[PPP] = prim[RHO]*1e-5;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;
}
