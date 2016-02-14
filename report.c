
#include "paul.h"

void report( struct domain * theDomain , double t ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Ng = theDomain->Ng;
   int * Nr = theDomain->Nr;
   int rank = theDomain->rank;
   int size = theDomain->size;
   MPI_Comm grid_comm = theDomain->theComm;

   int jmin = Ng;
   int jmax = Nt-Ng;
   if( rank==0 ) jmin = 0;
   if( rank==size-1 ) jmax = Nt;

   double L1 = 0.0;
   double uSum = 0.0;
   double ESum = 0.0;
   double EJet = 0.0;
   double XSum = 0.0;
   double MSum = 0.0;
   double Ni   = 0.0;
   double uMax = 0.0;
   double gh_max = 0.0;
   double rMax = 0.0;
   double rMin = HUGE_VAL;
   int i,j;
   for( j=jmin ; j<jmax ; ++j ){
      int jk = j;
      for( i=0 ; i<Nr[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         double ur = c->prim[UU1];
         double up = c->prim[UU2];
         double X = 0.0;
         if( NUM_Q > NUM_C ) X = c->prim[NUM_C];
         double u = sqrt(ur*ur+up*up);
         double E = c->cons[TAU];
         double M = c->cons[DEN];
         double h = 1.+4.*c->prim[PPP]/c->prim[RHO];
         double gam = sqrt( 1. + u*u );
         uSum += u*E;
         ESum += E;
         double u_eff = sqrt( gam*h*gam*h - 1.0 );
         if( u_eff > 1.0 ) EJet += E;
         XSum += X*E;
         MSum += M;
         Ni += X*M;
         if( uMax < u ) uMax = u;
         if( gh_max < gam*h ) gh_max = gam*h;
      }
   }

   MPI_Allreduce( MPI_IN_PLACE , &uMax   , 1 , MPI_DOUBLE , MPI_MAX , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &gh_max , 1 , MPI_DOUBLE , MPI_MAX , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &uSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &ESum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &EJet , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &XSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &MSum , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Ni   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   double uAv = uSum/ESum;

   for( j=jmin ; j<jmax ; ++j ){
      int jk = j;
      for( i=0 ; i<Nr[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         double ur = c->prim[UU1];
         double up = c->prim[UU2];
         double u = sqrt(ur*ur+up*up);
         double r = c->riph;
         if( rMax < r && u>.5*uAv ) rMax = r;
         if( rMin > r && u>.5*uAv ) rMin = r;
         L1 += fabs(c->prim[PPP]/pow(c->prim[RHO],5./3.)-1.)*c->dr;
      }
   }

   MPI_Allreduce( MPI_IN_PLACE , &rMax , 1 , MPI_DOUBLE , MPI_MAX , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &rMin , 1 , MPI_DOUBLE , MPI_MIN , grid_comm );

   double v_phot = 0.0;
   if( rank==0 ){
      int jk = 0;
      double tau = 0.0;
      for( i=Nr[jk]-1 ; i>=0 && tau<1. ; --i ){
         struct cell * c = &(theCells[jk][i]);
         double rho = c->prim[RHO];
         double kappa = 3.*t*t;
         double dr = c->dr;
         double r = c->riph;
         tau += rho*kappa*dr;
         double ur = c->prim[UU1];
         double up = c->prim[UU2];
         double gam = sqrt(1.+ur*ur+up*up);
         v_phot = ur/gam;
      }    
   }

   if( rank==0 ){
      FILE * rFile = fopen("report.dat","a");
      fprintf(rFile,"%e %e %e %e %e %e %e %e %e %e %e\n",t,rMax,rMin,uAv,uMax,ESum,MSum,EJet,gh_max,v_phot,Ni);
      fclose(rFile);
   }

}
