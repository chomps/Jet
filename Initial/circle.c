
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   double z0 = 0.6;
   double x0 = 0.0;
   double r  = x[0];
   double th = x[1];
   double ph = x[2];
   double x1 = r*sin(th)*cos(ph)-x0;
   double y1 = r*sin(th)*sin(ph);
   double z1 = r*cos(th)-z0;
   double rr = sqrt(x1*x1+y1*y1+z1*z1);
   double rho,Pp;
   if( rr < .25 ){
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
