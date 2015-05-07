
#include "../paul.h"

static double R_MIN = 0.0; 
static double GAMMA_LAW = 0.0; 

void setICparams( struct domain * theDomain ){
   R_MIN = theDomain->theParList.rmin; 
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index; 
}

void initial( double * prim , double * x ){

   double r  = x[0];
   //double th = x[1];

   double gam = GAMMA_LAW;

   double R = 1.0;
   double rho0 = 1.0;
   double rhoc = .5*rho0;
   double Pmin = 1e-6;
   //double rhomin = 1e-5*rho0;

   double rho = rho0*( 1. - r*r/R/R );
   if( rho < 0.0 ) rho = 0.0;

   double Mass = (0.115)*( 4.*M_PI*R*R*R*rho0 );
   double Msol = Mass/4.77;

//   double Bethe = .000559*Msol; //10^51 erg / Msol c^2


   double rexp   = 5.*R_MIN; 
   double E = .0033*Mass;//200.*Bethe;//.0033*Mass;
   double e0 = E/pow(sqrt(M_PI)*rexp,3.);
   double ex = e0*exp(-r*r/rexp/rexp);
   double Pexp = (gam-1.)*ex;

   double Msol_py = 1.49e-8*Msol/R;
   double Mdot = 1e-5*Msol_py;
   double vwind = .012;

   double wind = Mdot/4./M_PI/vwind/r/r;


   if( rho < rhoc ){
      double XX = rho/rhoc - .5;
      if( XX < 0.0 ) XX = 0.0;
      rho = rhoc*sqrt(2.*XX);
   }

   rho += wind;
   //if( rho < rhomin ) rho = rhomin;
   double Pp = Pexp + Pmin*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}
