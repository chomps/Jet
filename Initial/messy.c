
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r = x[0];

   double r0 = 1.0;//1.3941e16;
   double rc = 0.1*r0;
   double t0 = 1.0;//86400*45.;
   double day = t0/45.;
   double rho_csm = 1.0;//1.3177e-17;

   double t_ej = 5*day;
   double rt = 0.31654*rc;
   double rho_ej  = 3.0307*rho_csm;

   double v0 = r0/t0;

   double Poverrho = 1e-5*v0*v0;//8.2e7;

   double v   = 0.0;
   double X   = 0.0;
   double rho = rho_csm;

   if( r < rc ){
      rho = rho_ej*pow(rc/r,10.);
      v = r/t_ej;
      X = 1.0;
      if( r < rt ) rho = rho_ej*pow(rc/rt,10.)*rt/r;
   }else if( r>2.*rc ){
      rho *= 1e-4;
   }
 
   prim[RHO] = rho;
   prim[PPP] = rho*Poverrho;
   prim[UU1] = v;
   prim[UU2] = 0.0;
   if( NUM_N > 0 ) prim[NUM_C] = X;

}
