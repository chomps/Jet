
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../paul.h"
 
void prim2cons( double * , double * , double , double , double );
void flux( double * , double * , double , double , double * );
void getUstar( double * , double * , double , double , double , double , double * );
void vel( double * , double * , double * , double * , double * , double * , double , double );
 
void riemann_r( struct cell * cL , struct cell * cR , double r , double th , double dAdt ){

   double primL[NUM_Q];
   double primR[NUM_Q];

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + .5*cL->gradr[q]*cL->dr;
      primR[q] = cR->prim[q] - .5*cR->gradr[q]*cR->dr;
   }

   double Sl,Sr,Ss;
   double n[3];
   n[0] = 1.0;   n[1] = 0.0;   n[2] = 0.0;

   vel( primL , primR , &Sl , &Sr , &Ss , n , r , th );

   double Fk[NUM_Q];
   double Uk[NUM_Q];

   double Flux[NUM_Q];

   if( MOVE_CELLS == C_WRIEMANN ) cL->wiph += Ss;
   double w = cL->wiph;

   if( w < Sl ){
      flux( primL , Fk , r , th , n );
      prim2cons( primL , Uk , r , th , 1.0 );
      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fk[q] - w*Uk[q];
      }
   }else if( w > Sr ){
      flux( primR , Fk , r , th , n );
      prim2cons( primR , Uk , r , th , 1.0 );
      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fk[q] - w*Uk[q];
      }
   }else{
      double Ustar[NUM_Q];
      if( w < Ss ){
         prim2cons( primL , Uk , r , th , 1.0 );
         getUstar( primL , Ustar , r , th , Sl , Ss , n );
         flux( primL , Fk , r , th , n );

         for( q=0 ; q<NUM_Q ; ++q ){
            Flux[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] ) - w*Ustar[q];
         }
      }else{
         prim2cons( primR , Uk , r , th , 1.0 );
         getUstar( primR , Ustar , r , th , Sr , Ss , n );
         flux( primR , Fk , r , th , n );

         for( q=0 ; q<NUM_Q ; ++q ){
            Flux[q] = Fk[q] + Sr*( Ustar[q] - Uk[q] ) - w*Ustar[q];
         }
      }
   }

   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dAdt;
      cR->cons[q] += Flux[q]*dAdt;
   }

}

void riemann_trans( struct face * F , double dt , int dim ){

   struct cell * cL = F->L;
   struct cell * cR = F->R;
   double dA        = F->dA;
   double dxL       = F->dxL;
   double dxR       = F->dxR;
   double r         = F->cm[0];
   double th        = F->cm[1];

   double primL[NUM_Q];
   double primR[NUM_Q];

   double rL = cL->riph - .5*cL->dr;
   double rR = cR->riph - .5*cR->dr;
   double drL = r - rL;
   double drR = rR - r;

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + cL->grad[q]*dxL + cL->gradr[q]*drL;
      primR[q] = cR->prim[q] - cR->grad[q]*dxR - cR->gradr[q]*drR;
   }
   

   double Sl,Sr,Ss;
   double n[3] = {0.0,0.0,0.0};
   n[dim] = 1.0;
 
   vel( primL , primR , &Sl , &Sr , &Ss , n , r , th );

   double Fk[NUM_Q];
   double Uk[NUM_Q];
   
   double Flux[NUM_Q];

   if( 0. < Sl ){
      flux( primL , Flux , r , th , n );
   }else if( 0. > Sr ){
      flux( primR , Flux , r , th , n );
   }else{
      double Ustar[NUM_Q];
      if( 0. < Ss ){
         prim2cons( primL , Uk , r , th , 1.0 );
         getUstar( primL , Ustar , r , th , Sl , Ss , n );
         flux( primL , Fk , r , th , n );

         for( q=0 ; q<NUM_Q ; ++q ){
            Flux[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] );
         }
      }else{
         prim2cons( primR , Uk , r , th , 1.0 );
         getUstar( primR , Ustar , r , th , Sr , Ss , n );
         flux( primR , Fk , r , th , n );

         for( q=0 ; q<NUM_Q ; ++q ){
            Flux[q] = Fk[q] + Sr*( Ustar[q] - Uk[q] );
         }
      }
   }

   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= dt*dA*Flux[q];
      cR->cons[q] += dt*dA*Flux[q];
   }
}

