
#include "../paul.h"
#include <string.h>

void initial( double * , double * );
double get_dV( double * , double * );
void cons2prim( double * , double * , double , double , double );

void boundary_r( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int j,k,q;
   int ABSORB_R0 = theDomain->theParList.Absorb_BC;
   
   for( j=0 ; j<Nt ; ++j ){
      double theta = .5*(t_jph[j]+t_jph[j-1]);
      for( k=0 ; k<Np ; ++k ){
         double phi = .5*(p_kph[k]+p_kph[k-1]);
         int jk = j+Nt*k;
         struct cell * c2 = &(theCells[jk][Nr[jk]-1]);
         double r = c2->riph - c2->dr;
         double x[3] = {r,theta,phi};
         initial( c2->prim , x );
         if( ABSORB_R0 ){
            struct cell * c3 = &(theCells[jk][1]);
            struct cell * c4 = &(theCells[jk][0]);
            for( q=0 ; q<NUM_Q ; ++q ){
               c4->prim[q] = c3->prim[q];
               if( c4->prim[UU1] > 0.0 ) c4->prim[UU1] = 0.0;
            }
         }
      }    
   }

}

double get_dV( double * , double * );
void prim2cons( double * , double * , double , double );

void boundary_trans( struct domain * theDomain , struct face * theFaces , int * nn , int dim ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   int copy_left  = 0;
   int copy_right = 0;

   int LRmin = 0; 
   int LRmax = 2;
   if( !copy_left  ) LRmin = 1;
   if( !copy_right ) LRmax = 1;

   if( dim==1 && (dim_rank[0]==0 || dim_rank[0]==dim_size[0]-1) ){
      int k,LR;
      for( LR=LRmin ; LR<LRmax ; ++LR ){
         if( ( dim_rank[0]==0 && LR==0 ) || ( dim_rank[0]==dim_size[0]-1 && LR==1 ) ){
            int j0 = 0;
            int j1 = 1;
            if( LR==1 ){
               j0 = Nt-1;
               j1 = Nt-2;
            }
            for( k=0 ; k<Np ; ++k ){
               int jk0 = j0+Nt*k;
               int jk1 = j1+Nt*k;
               Nr[jk0] = Nr[jk1];
               theCells[jk0] = realloc( theCells[jk0] , Nr[jk0]*sizeof(struct cell) );
               memcpy( &(theCells[jk0][0]) , &(theCells[jk1][0]) , Nr[jk0]*sizeof(struct cell) );
               int i;
               for( i=0 ; i<Nr[jk0] ; ++i ){
                  struct cell * c = &( theCells[jk0][i] );
                  c->prim[UU2] *= .5;
                  double rp = c->riph;
                  double rm = 0.0;
                  if( i>0 ) rm = theCells[jk0][i-1].riph;
                  //double r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
                  double xp[3] = {rp,t_jph[j0]  ,p_kph[k]  };
                  double xm[3] = {rm,t_jph[j0-1],p_kph[k-1]};
                  double dV = get_dV( xp , xm );
                  xp[1] = t_jph[j1  ];
                  xm[1] = t_jph[j1-1];
                  double dV2 = get_dV( xp , xm );
                  //prim2cons( c->prim , c->cons , r , dV );
                  int q;
                  for( q=0 ; q<NUM_Q ; ++q ){
                     c->cons[q]   *= dV/dV2;
                     c->RKcons[q] *= dV/dV2;
                  }
               }
            }  
         } 
      }
   }

}

