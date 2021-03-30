
#include "../paul.h"

static double T_MIN = 0.0;

void setICparams( struct domain * theDomain ){
   T_MIN = theDomain->theParList.t_min;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double rho0 = 1.0;
   double Pmin = 1e-5;

   double k = 2.0;

   double E = 1.0;
   double t = T_MIN;
   double l = pow(E/rho0,1./3.);
   double rho_out = rho0*pow(t/l,-k);
   double G = sqrt((17.-4.*k)/8./M_PI)*pow(E/rho_out/t/t/t,.5);
   double m = 3.-k;

//   double R = ( 1. - 1./(2.*(m+1.)*G*G) )*t;

   double chi = (1.+2.*(m+1.)*G*G)*(1.-r/t);

   double f = pow(chi,-(17.-4.*k)/(12.-3.*k));
   double g = 1./chi;
   double h = pow(chi,-(7.-2.*k)/(4.-k));

   double Pp  = (2./3.)*rho_out*G*G*f;//pow(chi,-17./12.);
   double gam = G*sqrt(.5*g); //sqrt(2.*chi);
   double rho = 2.*rho_out*G*G*h/gam; //pow(chi,-7./4.)/gam;
   
   if( chi < 1.0 ){//r > R ){
      rho = rho0*pow(r/l,-k);
      Pp = Pmin;
      gam = 0.0;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = gam;
   prim[UU2] = 0.0;
   //printf("gam = %e chi = %e G/r2 = %e R = %e\n",gam,chi,G/sqrt(2.),R);
}
