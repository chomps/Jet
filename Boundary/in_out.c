
#include "../paul.h"

double get_dV( double * , double * ); 
void prim2cons( double * , double * , double , double , double );
void initial( double * , double * );

void boundary_r( struct domain * theDomain ){

   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   struct cell ** theCells = theDomain->theCells;

   double t = theDomain->t;
   double vw = 1.0;

   double Amp = theDomain->theParList.Explosion_Energy;
   double wavenumber = theDomain->theParList.Gam_0;
 
   int j,k;
   for( j=0 ; j<Nt ; ++j ){
      double tp = t_jph[j];
      double tm = t_jph[j-1];
      for( k=0 ; k<Np ; ++k ){
         double pp = p_kph[k];
         double pm = p_kph[k-1];
         int jk = j+Nt*k;
         struct cell * c1 = &(theCells[jk][0]);
         struct cell * c2    = &(theCells[jk][Nr[jk]-1]);
         struct cell * c2_c  = &(theCells[jk][Nr[jk]-2]);

         double r1 = c1->riph;
         double x[3] = {r1,.5*(tp+tm),.5*(pp+pm)};
         double xp[3] = {r1, tp , pp};
         double xm[3] = {0.0,tm,pm};
         double dV = get_dV( xp , xm );
         initial( c1->prim , x );

double phase1 = .5*(tp+tm);
double phase2 = vw*t/r1;
double k_p = wavenumber;
c1->prim[RHO] *= 1. + Amp*sin(2.5*k_p*phase1)*sin(k_p*phase2);
c1->prim[RHO] *= 1. + .99*2.*((double)rand()/(double)RAND_MAX-.5);
 
         prim2cons( c1->prim , c1->cons , r1 , .5*(tp+tm) , dV );
         int q;
         for( q=0 ; q<NUM_Q ; ++q ){
            c1->RKcons[q] = c1->cons[q];
            c2->prim[q] = c2_c->prim[q];
         }
      }    
   }
}


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
   int copy_right = 1; 

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


