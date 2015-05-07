
#include "../paul.h"

static double T_MIN = 0.0;

void setICparams( struct domain * theDomain ){
   T_MIN = theDomain->theParList.t_min;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double rho0 = 1.0;
   double Pmin = 1e-5;


   double E = 1.0;
   double t = T_MIN;
   double l = pow(E/rho0,1./3.);
   double G = sqrt(17./8./M_PI)*pow(t/l,-1.5);

/*
   double G = 100.0;
   double t = 1.0;
*/
   double R = (1.-1./8./G/G)*t;

   double chi = (1.+8.*G*G)*(1.-r/t);

   double Pp  = (2./3.)*rho0*G*G*pow(chi,-17./12.);
   double gam = G/sqrt(2.*chi);
   double rho = 2.*rho0*G*G*pow(chi,-7./4.)/gam;
   
   if( chi < 1.0 ){//r > R ){
      rho = rho0;
      Pp = Pmin;
      gam = 0.0;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = gam;
   prim[UU2] = 0.0;
   //printf("gam = %e chi = %e G/r2 = %e R = %e\n",gam,chi,G/sqrt(2.),R);
}
