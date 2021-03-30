
#include "../paul.h"

static double T_MIN     = 0.0;

void setICparams( struct domain * theDomain ){
    T_MIN    = theDomain->theParList.t_min;
}

double f( double r , double a , double b , double M ){
   return( (M/M_PI/a)/( (r-b)*(r-b) + 1./a/a ) );
}

void initial( double * prim , double * x ){

   double w = 1.0;

   double v_max = 0.3; //0.2;

   double M0 = 1.0;
   double M_ext = 1e-10*5e-4/.04*M0;
   double R_ext = 1.0;
   double t = T_MIN;
   double R0 = v_max*t;
   double r  = x[0];
   double th = x[1];
   
   r *= pow(w,2./3.)*sqrt( cos(th)*cos(th) + sin(th)*sin(th)/w/w );

   double R1 = .06*R0;
   double rhoc = M0/pow(R0,3.)/.0102913;  ///.003479;

//	2	.0221136
//	2.5	.0102913
//	3	.00554687
//	3.5	.00347901

   double rhostar = rhoc/( 1. + pow(r/R1,2.5) )*pow(1.-r/R0,.75);
   if( r>R0 ) rhostar = 0.0;;

   double rho_trans = 1e-8;
   //double rho_ISM = 1e-27;
   double rho_ISM = 1e-20;
   double r_ext = x[0];

//   double rho_ext = M_ext/4./M_PI/r_ext/r_ext/R_ext*exp(-r_ext/R_ext);
   double rho_ext = rho_trans*exp(-r_ext/R0);
//   double rho_ext = M_ext/4./M_PI/r_ext/r_ext/R_ext/(1. + pow(r_ext/R_ext,8.));

//   double rho_ext_R = M_ext/4./M_PI/pow(R_ext,3.);
//   double rho_trans = 0.0;
//   double rho_trans = sqrt(rho_ext_R*rho_ISM);
//   rho_trans = sqrt( rho_trans*rho_ext_R );

//   rho_ISM += rho_trans/(1.+pow(r_ext/R_ext,3.) );
   rho_ext += rho_ISM;
   double rho = rhostar + rho_ext;

   double v = v_max*(x[0]/R0)*rhostar/rho;
   double u = v/sqrt(1.-v*v);

   double Pp = 1e-6*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = u;
   prim[UU2] = 0.0;

   double X = rho_ISM/rho;
   if( X<.5 ) X=0.; else X = 1.0;

   if( NUM_N > 0 ) prim[NUM_C] = X;

}
