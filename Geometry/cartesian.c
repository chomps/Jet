
#include "../paul.h"

double get_dL( double * xp , double * xm , int dim ){
   return( xp[dim]-xm[dim] );
}

double get_dA( double * xp , double * xm , int dim ){
   double dr  = xp[0]-xm[0];
   double dth = xp[1]-xm[1];
   double dph = xp[2]-xm[2];
   if( dim==0 ) return( dth*dph );
   else if( dim==1 ) return( dr*dph );
   else return( dr*dth );
}

double get_dV( double * xp , double * xm ){
   double dr  = xp[0]-xm[0];
   double dth = xp[1]-xm[1];
   double dph = xp[2]-xm[2];
   return( dr*dth*dph );
}
