
#include "../paul.h"

void initial( double * , double * );

void boundary_r( struct cell ** theCells ){
   int j,k;
   for( j=0 ; j<Nt ; ++j ){
      double tp = t_jph[j];
      double tm = t_jph[j-1];
      for( k=0 ; k<Np ; ++k ){
         double pp = p_kph[k];
         double pm = p_kph[k-1];
         int jk = j+Nt*k;
         struct cell * c1 = &(theCells[jk][0]);
         struct cell * c2 = &(theCells[jk][Nr[jk]-1]);
         double r1 = .5*c1->riph;
         double r2 = c2->riph-.5*c2->dr;
         double x[3] = {r1,.5*(tp+tm),.5*(pp+pm)};
         initial( c1->prim , x );
         x[0] = r2;
         initial( c2->prim , x );
      }    
   }
}

void boundary_trans( struct cell ** theCells , struct face * theFaces , int * nn , int dim ){
   if( dim==1 ){
      int lmost = dim_rank[0]==0;
      int rmost = dim_rank[0]==dim_size[0]-1;
      int i,k,LR;
      for( LR=0 ; LR<2 ; ++LR ){
         if( (LR==0 && lmost) || (LR==1 && rmost) ){
            int j=0;
            if( LR==1 ) j=Nt-1;
            double th = .5*(t_jph[j]+t_jph[j-1]);
            for( k=0 ; k<Np ; ++k ){
               int jk = j+Nt*k;
               for( i=0 ; i<Nr[jk] ; ++i ){
                  struct cell * c = &(theCells[jk][i]);
                  double r  = c->riph-.5*c->dr;
                  double ph = .5*(p_kph[k]+p_kph[k-1]);
                  double x[3] = {r,th,ph};
                  initial( c->prim , x );
               }
            }
         }
      }
   }
}

