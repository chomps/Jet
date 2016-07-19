
#include "../paul.h"

static double delt0 = 0.0;

void setICparams( struct domain * theDomain ){
   srand(theDomain->rank);
   delt0 = theDomain->theParList.Explosion_Energy;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double R      = 1.0;  //300 pc
   double rho0   = 1.0;
   double vw     = 1.0;
   double gam = 5./3.;
   double P0 = rho0*vw*vw/gam;

   double mdot = rho0*vw*R*R;
   double s    = P0/pow(rho0,gam);
   double k    = .5*vw*vw + P0/rho0/(gam-1.);

//   mdot *= 1. + 1.5*((double)rand()/(double)RAND_MAX-.5);
 
   double v = vw;

   int j;
   for( j=0 ; j<30; ++j ){
      double v2 = fabs( 2.*k - 2.*s/(gam-1.)*pow(mdot/v/r/r,gam-1.) );
      v = sqrt(v2);
   }

   double rho = mdot/v/r/r;
   double Pp = s*pow(rho,gam);

   //rho *= 1. + 1.5*((double)rand()/(double)RAND_MAX-.5);
   double rho_adjust = sinh(delt0/2.)/(delt0/2.);
   double dx = delt0*((double)rand()/(double)RAND_MAX - .5);
   rho *= exp(dx)/rho_adjust;
   //rho *= 1. + 10.*delt0*((double)rand()/(double)RAND_MAX-.5);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = v;
   prim[UU2] = 0.0;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ) prim[q] = 0.0;
}


