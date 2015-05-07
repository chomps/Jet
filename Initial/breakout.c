
#include "../paul.h"

static double R_MIN = 0.0;
static double GAMMA_LAW = 0.0;
static double EXP_E = 0.0;

void setICparams( struct domain * theDomain ){
   R_MIN = theDomain->theParList.rmin; 
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index; 
   EXP_E = theDomain->theParList.Explosion_Energy;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double rho,Pp;

   double aspect = 1.0;

//   double rr = r;
   double stretch = sqrt( 1. + (1./aspect/aspect-1.)*sin(th)*sin(th));
   stretch *= pow(aspect,2./3.);
   r *= stretch;
 
   double R = 1.0;//sqrt(1.+(1./aspect/aspect-1.)*sin(th)*sin(th));
   double n = 3.85;
   double k = 1.9;

   double rc     = .5*R;

   double rexp   = 5.*R_MIN; 
   double rho0 = 1.0; 
   double Mass = 1.17*rho0*(4.*M_PI/(3.-k))*pow(.5*R,3.);
   double E =  EXP_E*Mass;
   double e0 = E/pow(sqrt(M_PI)*rexp,3.);
   double ex = e0*exp(-r*r/rexp/rexp);
   double Pexp = (GAMMA_LAW-1.)*ex;

   double Msol = Mass/4.77;

   double Msol_py = 1.49e-8*Msol/R;
   double Mdot = 1e-5*Msol_py;
   double vwind = .012;

   double wind = Mdot/4./M_PI/vwind/r/r;

   if( r < rc ){
      rho = rho0*pow(rc/r,k);
   }else if( r < R ){ 
      rho = rho0*(pow((R/r-1.0)/(R/rc-1.0),n) + wind);
   }else{
      rho = rho0*wind;
   }    

   Pp = Pexp + 1e-6*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}
