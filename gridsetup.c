
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "paul.h"

int getN0( int drank , int dsize , int dnum ){
   int N0 = (dnum*drank)/dsize;
   return(N0);
}

void gridSetup( struct domain * theDomain ){

   int Ng = NUM_G;
   theDomain->Ng = Ng;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int NUM_T = theDomain->theParList.Num_T;
   int NUM_P = theDomain->theParList.Num_P;

   double THETA_MIN = theDomain->theParList.thmin;
   double THETA_MAX = theDomain->theParList.thmax;
   double PHI_MAX   = theDomain->theParList.phimax;

   int N0t = getN0( dim_rank[0] , dim_size[0] , NUM_T );
   int N1t = getN0( dim_rank[0]+1 , dim_size[0] , NUM_T );
   int Nt = N1t-N0t;

   int N0p = getN0( dim_rank[1] , dim_size[1] , NUM_P );
   int N1p = getN0( dim_rank[1]+1 , dim_size[1] , NUM_P );
   int Np = N1p-N0p;

   if( dim_rank[0] != 0 ){ Nt += Ng; N0t -= Ng; }
   if( dim_rank[0] != dim_size[0]-1 ) Nt += Ng;
   if( NUM_P != 1 ) Np += 2*Ng;

   theDomain->Nt = Nt;
   theDomain->Np = Np;

//printf("rank %d:  N0t = %d, N1t = %d, Nt = %d.\n",theDomain->rank,N0t,N1t,Nt);

   theDomain->Nr    = (int *)    malloc( Np*Nt*sizeof(int) );
   theDomain->t_jph = (double *) malloc( (Nt+1)*sizeof(double) );
   theDomain->p_kph = (double *) malloc( (Np+1)*sizeof(double) );

   ++(theDomain->t_jph);
   ++(theDomain->p_kph);

   int j,k;

//   double dth = (THETA_MAX-THETA_MIN)/(double)NUM_T;
//   double t0 = THETA_MIN + (double)N0t*dth;
   double dx = 1./NUM_T;
   double x0 = (double)N0t*dx;
   double th0 = .4;
   double thmax = THETA_MAX;
   double k0 = log( tan(.5*(th0+thmax))/tan(.5*th0) );

   for( j=-1 ; j<Nt ; ++j ){
      double x = x0 + ((double)j+1.)*dx;
//      double theta = THETA_MIN + 2.*atan( tan(.5*th0)*exp(k0*x) )-th0;//(THETA_MAX-THETA_MIN)*x;
      double theta = THETA_MIN + (THETA_MAX-THETA_MIN)*x;
      theDomain->t_jph[j] = theta;//t0 + ((double)j+1.)*dth;
   }
   double dphi = PHI_MAX/(double)NUM_P;
   double p0 = (double)N0p*dphi;
   for( k=-1 ; k<Np ; ++k ){
      theDomain->p_kph[k] = p0 + ((double)k+1.)*dphi;
   }

   for( j=0 ; j<Nt ; ++j ){
      double dth_j = 1.001*(theDomain->t_jph[j]-theDomain->t_jph[j-1]);
      double dth = (THETA_MAX-THETA_MIN)/(double)NUM_T;
      double nrj = theDomain->theParList.Num_R*dth/dth_j;
      for( k=0 ; k<Np ; ++k ){
         theDomain->Nr[j+k*Nt] = (int)nrj;
      }
   }

}


