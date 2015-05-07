
#include "../paul.h"

double get_dL( double * xp , double * xm , double dx , int dim ){
   double r  = .5*(xp[0]+xm[0]);
   double th = .5*(xp[1]+xm[1]);
   if( dim==0 ) return( xp[0]-xm[0] );
   else if( dim==1 ) return( r*(xp[1]-xm[1]) );
   else return( r*sin(th)*( xp[2]-xm[2] ) );
}

double get_dA( double * xp , double * xm , int dim ){
   double r  = .5*(xp[0]+xm[0]);
   double th = .5*(xp[1]+xm[1]);
   double dr  = xp[0]-xm[0];
   double dth = xp[1]-xm[1];
   double dph = xp[2]-xm[2];
   double sinth = sin(th)*sin(.5*dth )/(.5*dth);
   if( dim==0 ) return( r*r*sinth*dth*dph );
   else if( dim==1 ) return( r*sin(th)*dr*dph );
   else return( r*dr*dth );
}

double get_dV( double * xp , double * xm ){
//   double r  = .5*(xp[0]+xm[0]);
   double th = .5*(xp[1]+xm[1]);
   double dr  = xp[0]-xm[0];
   double dth = xp[1]-xm[1];
   double dph = xp[2]-xm[2];

   double r2    = (xp[0]*xp[0]+xm[0]*xm[0]+xp[0]*xm[0])/3.;
   double sinth = sin(th)*( sin(.5*dth)/(.5*dth) );

   return( r2*sinth*dr*dth*dph );

}
