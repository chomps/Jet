
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   //double x0 = 0.75;
   //double R  = 0.4;
   double A  = 3.0;//7e2;

   double r  = x[0];
   //double th = x[1];
   //double x1 = r*sin(th);
   //double y1 = r*cos(th);//-x0;
   //double rr = sqrt(x1*x1+y1*y1);
   double rho,Pp,delta;
   delta = 0.0;
   //if( rr < R ){
      delta = A*exp(-80.*r*r);//exp(1./(rr*rr-R*R));
   //}
   rho = 1.0*(1.+delta);;
   Pp = pow(rho,5./3.);
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}

