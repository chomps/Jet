
#include "paul.h"

//static double r0;

double get_dV( double * , double * );

void setCoolingParams( struct domain * theDomain ){

//   r0   = theDomain->theParList.Nozzle_r0;

}

void cool_src( double * prim , double * cons , double dVdt ){

   double T0 = 3e-3; //1e5 K
   double Tmin = 0.05*T0;

   double Q = 0.0;
   double Pp  = prim[PPP];
   double rho = prim[RHO];

   double T = Pp/rho;
   if( T > Tmin ){
      Q = -5.0*rho*rho*(T-Tmin)/T0*pow(T/T0,-1.7);
   }

   cons[TAU] += Q*dVdt;

}

void add_cooling( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double t = theDomain->t;

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
            double xp[3] = {rp,thp,php};
            double xm[3] = {rm,thm,phm};
            double dV = get_dV(xp,xm);
            cool_src( c->prim , c->cons , dV*dt );
         }
      }
   }

}

