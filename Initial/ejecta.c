
#include "../paul.h"

static double t = 0.0;

void setICparams( struct domain * theDomain ){
   t = theDomain->theParList.t_min;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double E_ej = 1.0;
   double M_ej = 1.0;

   double n = 10.0;
   double d = 1.0;

//Ejecta Density
   double vt2 = E_ej/M_ej*2.*(5.-d)*(n-5.)/((3.-d)*(n-3.));
   double vt = sqrt(vt2);
   double R = vt*t;

   double zeta = ( (n-3.)*(3.-d)/(n-d) )/4./M_PI;
   double rho_T = zeta*M_ej/pow(R,3.);

   double rho_ej = rho_T*pow(r/R,-d);
   if( r>R ) rho_ej = rho_T*pow(r/R,-n);

   double rho = rho_ej;
   double Pp = 1e-10*rho*vt2;
   double v = r/t;
   if( r>2.*R ) v *= exp(-pow((r-2.*R)/2./R,4.));

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = v;
   prim[UU2] = 0.0;

}


