
#include "../paul.h"

static double T_MIN = 0.0; 

void setICparams( struct domain * theDomain ){
   T_MIN = theDomain->theParList.t_min;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double k = 0.0;

   double t   = T_MIN;
   double th_j = M_PI;
   double Gam = 30.0;
   double delta = 1.0;

   double E   = 1.0;
   double rho0 = 1.0;
   double r_Sedov = pow(E/rho0,1./3.);

   double K = 1e-8;
   double gmax = Gam*delta/(delta-0.5);

   double R = t*sqrt(1.- 1./gmax/gmax);
   double rhomax = delta/2./M_PI*E/pow(t,3.);

   double rho = rhomax*pow( (1.-R/t)/(1.-r/t) , delta );
   double v = r/t;

   double X = 1.0;

   if( r>R || th>th_j ){
      X = 0.0;
   }

   rho = X*rho + (1.-X)*( rho0 * pow(r_Sedov/r,k) );
   double w = X*exp(-pow(2.*(r-R)/R,4.));
   v   = w*v;

   double u = v/sqrt(1.-v*v);
   double Pp = K*pow(rho,4./3.);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = u;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = X;

}
