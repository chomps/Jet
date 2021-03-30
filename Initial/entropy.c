
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double R0 = 5.0;
   double t0 = 0.0;//1.45;
   double p0 = 0.0;//0.15;
   
   double x0 = R0*sin(t0)*cos(p0);
   double y0 = R0*sin(t0)*sin(p0);
   double z0 = R0*cos(t0);

   double A  = 30.0;//7e2;

   double r  = x[0];
   double th = x[1];
   double ph = x[2];
   double x1 = r*sin(th)*cos(ph)-x0;
   double y1 = r*sin(th)*sin(ph)-y0;
   double z1 = r*cos(th)-z0;
   double rr = sqrt(x1*x1+y1*y1+z1*z1);
   double rho,Pp,delta;
   delta = 0.0;
   //if( rr < R ){
      delta = A*exp(-3.*rr*rr);//exp(1./(rr*rr-R*R));
   //}
   rho = 1.0*(1.+delta);;
   Pp = pow(rho,4./3.);
   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;
   if(NUM_C > UU3) prim[UU3] = 0.0;

}

