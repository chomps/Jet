
#include "../paul.h"

static double R_MIN = 0.0;
static double t = 0.0;
static double hr = 0.0;

void setICparams( struct domain * theDomain ){
   t = theDomain->theParList.t_min;
   R_MIN = theDomain->theParList.rmin;
   hr = theDomain->theParList.Explosion_Energy;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double E_ej = 100.0;
   double M_ej = 100.0;
   double M_disk = 1.0;

   double A_wind = 1e-8;
   double R_disk = 1.0;

//Wind Density
   double rho_wind = A_wind/r/r;

//Ejecta Density
   double d = 1.0;
   double n = 10.0;
   double vt2 = E_ej/M_ej*2.*(5.-d)*(n-5.)/((3.-d)*(n-3.));
   double vt = sqrt(vt2);
   double zeta = ( (n-3.)*(3.-d)/(n-d) )/4./M_PI;
   double rho_T = zeta*M_ej/pow(vt*t,3.);

   double rho_ej = rho_T*pow(r/vt/t,-d);
   if( r>vt*t ) rho_ej = rho_T*pow(r/vt/t,-n);
   if( r>10.*vt*t ) rho_ej = 0.0;

//Disk Density
   
   double rho_disk = sqrt(1./hr/hr+.25*M_PI)*(M_disk/pow(R_disk,3.)*pow(r/R_disk,-2.))*exp(2./hr/hr*(sin(th)-1.));
   rho_disk *= exp(-pow(r/R_disk,4.))/10.2;  ///10.4;

   double rho = rho_disk + rho_wind + rho_ej;
   double Pp = 1e-3*rho;

   double v = r/t*(rho_ej/rho);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = v;
   prim[UU2] = 0.0;

}
