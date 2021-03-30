
#include "../paul.h"

static double delt0 = 0.0;

void setICparams( struct domain * theDomain ){
   srand(theDomain->rank);
   rand();
   delt0 = theDomain->theParList.Explosion_Energy;
}

double f_solve( double v , double k , double s , double mdot , double r , double gam ){
   return( -v*v + 2.*k - 2.*s/(gam-1.)*pow(mdot/v/r/r,gam-1.));
}

double dfdx_solve( double v , double s , double mdot , double r , double gam ){
   return( -2.*v + 2.*s*pow( mdot/v/r/r , gam-1. )/v );
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

   double v = vw;

   int j;
   for( j=0 ; j<4 ; ++j ){
      double f = f_solve( v , k , s , mdot , r , gam );
      double dfdx = dfdx_solve( v , s , mdot , r , gam );
      v -= f/dfdx;
   }

   double rho = mdot/v/r/r;
   double Pp = s*pow(rho,gam);

   //rho *= 1. + 1.5*((double)rand()/(double)RAND_MAX-.5);
//   double rho_adjust = sinh(delt0/2.)/(delt0/2.);
//   double dx = delt0*((double)rand()/(double)RAND_MAX - .5);
//   rho *= exp(dx)/rho_adjust;
   //rho *= 1. + 10.*delt0*((double)rand()/(double)RAND_MAX-.5);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = v;
   prim[UU2] = 0.0;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ) prim[q] = 0.0;
}


