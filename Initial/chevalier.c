
#include "../paul.h"

static double T_MIN = 0.0; 
static double delta = 0.0;

void setICparams( struct domain * theDomain ){
   T_MIN = theDomain->theParList.t_min;
   delta = theDomain->theParList.Explosion_Energy;
}
 
void initial( double * prim , double * x ){

   double n = 7.0;
   double s = 2.0; 

   double r  = x[0];
   double t   = T_MIN;

   double g = 1.0;
   double q = 1.0;

   double rho1 = pow(r/t/g,-n)*pow(t,-3.);
   double rho2 = q*pow(r,-s);

   double R0 = pow( pow(t,n-3.)*pow(g,n)/q , 1./(n-s) );

   double rho = rho1+rho2;
   double X   = rho1*rho1/(rho1*rho1+rho2*rho2);
   double v   = (r/t)*X;

   double r1 = 0.015*R0;
//   double r1 = 0.007*R0;
//   double r1 = 0.002*R0;
//   double r1 = 0.1*R0;
   if( r<r1 ){
      //rho = pow(r1/t/g,-n)*pow(t,-3.)*pow(r/r1,-1.);
      rho = pow(r1/t/g,-n)*pow(t,-3.)*pow(r/r1,-n*r/r1);
//      v *= pow(r/r1,20.);
   }

   if( X < 0.5 ){
      X = 0.0;
   }else{
      X = 1.0;
   }

   double k = 50.0;
   double th = x[1];
   double ptb = delta*sin(k*log(r))*sin(k*th);

   double v0 = pow( r , (s-3.)/(n-3.) );
   if( v0 > r/t ) v0 = r/t;
   double Mach = 30.; //1e3;
   double cs = v0/Mach;

   double Pmin = rho*cs*cs;

   prim[RHO] = rho*exp(ptb);
   prim[PPP] = Pmin;
   prim[UU1] = v;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = X;

}
