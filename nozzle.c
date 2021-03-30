
#include "paul.h"

//static double r0;
static double Pow;
static double th0;
static double eta0;
static double gam0;
static double tjet;
static double tmin;

void setNozzleParams( struct domain * theDomain ){

//   r0   = theDomain->theParList.Nozzle_r0;
   Pow  = theDomain->theParList.Nozzle_Power;
   th0  = theDomain->theParList.Nozzle_th0; 
   eta0 = theDomain->theParList.Nozzle_Eta;
   gam0 = theDomain->theParList.Nozzle_Gamma;
   tjet = theDomain->theParList.Nozzle_Time;
   tmin = theDomain->theParList.t_min;

}

void noz_src( double * cons , double dVdt , double r , double theta , double t , double r_min ){

   double r0,v,Vol,f,eta;
   int wind = 0;
   if( !wind ){

   //double r0 = 5.*r_min;
      r0 = 2.*r_min;
      v = sqrt(1.-1./gam0/gam0);

      Vol = pow( sqrt(2.*M_PI)*r0 , 3. )*( 1. - exp(-2./th0/th0) )*th0*th0/sqrt(.5*M_PI);
      f = (r/r0)*exp(-.5*r*r/r0/r0)*exp( (fabs(cos(theta))-1.) /th0/th0)/Vol;
      eta = eta0;
   //f *= exp(-(t-tmin)/tjet);
   }else{
   
      r0 = 2.*r_min;
      v = 0.1;
      Vol = 8.*M_PI*pow( r0 , 3. );
      f = (r/r0)*exp(-.5*r*r/r0/r0)/Vol;
      eta = .5*v*v;
   }
//      f *= pow(1.+(t-tmin)/tjet,-2.);


   double SE = Pow*f;
   double SS = v*SE;
   double SM = SE/eta;

   cons[DEN] += SM*dVdt;
   cons[SS1] += SS*dVdt;// *cos(theta);
   //cons[SS2] += SS*dVdt*sin(theta)*(-r);
   cons[TAU] += SE*dVdt;

   if( NUM_N > 0 ){
      cons[NUM_C+0] += SM*dVdt;
   }


//   int q;
//   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
//      cons[q] += SM*dVdt;
//   }




/*
   double r0 = 3.*r_min;
   double Vol = 4./3.*M_PI*pow(r0,3.);
   double Q = 0.0;
   double v_wind = 50.0;
   if( r<r0 ) Q = Pow/Vol*pow(1.+t/tjet,-2.);
   cons[TAU] += Q*dVdt;
   cons[DEN] += Q*dVdt/eta0;
//   cons[SS1] += 2.*Q/v_wind*dVdt;
//   cons[DEN] += 2.*Q/v_wind/v_wind*dVdt;
//   int q;
//   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
//      cons[q] += 2.*Q/v_wind/v_wind*dVdt;
//   }

*/

/*
   double Rmin = 2.47*7e10;
   double Rmax = 5e11; //0.05;
   double Vol = 4./3.*M_PI*( pow(Rmax,3.) - pow(Rmin,3.) );
   double T = 1e7; //0.1;
   //double Pow = 3e41/(1.+pow(t/3e5,1.5));
//4.5e40
   double Pow = 1e41;

   double Q = 0.0;
   if( r<Rmax && r>Rmin && t<T ) Q = Pow/Vol;

   cons[TAU] += Q*dVdt;
*/

/*
   double R = 0.006;
   Vol = 4./3.*M_PI*pow(R,3.);
   double T = 1e10; 
   double E = 1.0;
//   double Pow = E/T;
//   Pow = 1e-4; 

   double Q = 0.0; 
   if( r<R && t<T ) Q = Pow/Vol;
   double eps = 1e-1;

   cons[TAU] += Q*dVdt;
   cons[DEN] += eps*Q*dVdt;
*/
}

void noz_set( double * prim , double r , double theta , double t , double r_min ){

   double r0 = 3.*r_min;
   if( r<r0 ){   
      double rho,Pp,gv,X;

      if( theta < th0 ){

         double Area = 2.*M_PI*r0*r0*(1.-cos(th0));

         Pp  = Pow/4./gam0/gam0/Area*exp(-t/tjet);
         rho = 4.*gam0*Pp/eta0;
         gv  = gam0;
         X   = 1.0;

      }else{

         rho = prim[RHO];
         gv  = 0.0;
         Pp  = prim[PPP];
         X   = 0.0;

      }

      prim[RHO] = rho;
      prim[UU1] = gv;
      prim[UU2] = 0.0;
      prim[PPP] = Pp;
      int q;
      for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
         prim[q] = X;
      }
   }

}

double get_dV( double * , double * );

void nozzle( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double t = theDomain->t;
   double r_min = theCells[0][0].riph;

   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
      double thp = t_jph[j];
      double thm = t_jph[j-1];
      double th = .5*(thp+thm);
      for( k=0 ; k<Np ; ++k ){
         double php = p_kph[k];
         double phm = p_kph[k-1];
         int jk = j+Nt*k;
         for( i=0 ; i<Nr[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double rp = c->riph;
            double rm = rp - c->dr;
            double r  = rp - .5*c->dr;
            double xp[3] = {rp,thp,php};
            double xm[3] = {rm,thm,phm};
            double dV = get_dV(xp,xm);
            noz_src( c->cons , dV*dt , r , th , t , r_min );
            //noz_set( c->prim , r , th );
         }
      }
   }

}

