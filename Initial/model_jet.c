
#include "../paul.h"

static double T_MIN     = 0.0; 
static double GAM_EXP   = 4.5;
static double GAM_BOOST = 6.0;

void setICparams( struct domain * theDomain ){
   T_MIN    = theDomain->theParList.t_min;
   GAM_EXP   = theDomain->theParList.Gam_0;
   GAM_BOOST = theDomain->theParList.Gam_Boost;
}
 
double shell( double x , double z , double t , double Gam , double G_boost , int * in ){

   double r = sqrt( x*x + z*z );

   double E   = 2.0/Gam/G_boost; 
   double rho0 = 1.0; 
   double r_Sedov = pow(E/rho0,1./3.);

   double gmax = Gam;

   double R = t*Gam/sqrt(1.+Gam*Gam);
   double rhomax = 1./2./M_PI*E/pow(t,3.);

   double rho = rhomax*(1.-R/t)/(1.-r/t);
   *in = 1;

   if( r>R ){ rho = rho0 ; *in=0 ; }
   
   return( rho );

}

void boost( double gam , double z , double t , double * znew , double * tnew ){
   double v = sqrt( 1. - 1./gam/gam );
   *znew = gam*(z-v*t);
   *tnew = gam*(t-v*z);
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double t   = T_MIN;

   //double K = 1e-10;
   double G_boost = GAM_BOOST;
   double G_int   = GAM_EXP;
   double zz,tt;
   int inside;
   boost( G_boost , r*cos(th) , t , &zz , &tt );
   double rho = shell( r*sin(th) , zz , tt , G_int , G_boost , &inside );

   if( !inside ) rho = pow(r,-2.);

   double ur  = 0.0;
   double X   = 0.0;
   if( inside ){
      double vr = r/t;
//      if( vr < 0.5 ) vr = 0.0;
      ur = vr/sqrt(1.-vr*vr);
      X = 1.0;
   }
   //double Pp = K*pow(rho,4./3.);
   double Pp = 1e-5*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = ur;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = X;

}
