
#include "../paul.h"

void prim2cons( double * , double * , double , double , double );
void flux( double * , double * , double , double * );
void vel( double * , double * , double * , double * , double * , double * , double , double );

void riemann_r( struct cell * cL , struct cell * cR, double r , double dAdt ){

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

   vel( primL , primR , &Sl , &Sr , &Ss , n , r );

   double Fl[NUM_Q];
   double Fr[NUM_Q];
   double Ul[NUM_Q];
   double Ur[NUM_Q];

   double Flux[NUM_Q];

   if( MOVE_CELLS == C_WRIEMANN )  cL->wiph += Ss;
   double w = cL->wiph;

   if( w < Sl ){
      flux( primL , Fl , r , n );
      prim2cons( primL , Ul , r , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fl[q] - w*Ul[q];
      }
   }else if( w > Sr ){
      flux( primR , Fr , r , n );
      prim2cons( primR , Ur , r , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fr[q] - w*Ur[q];
      }
   }else{

      double Fstar;
      double Ustar;

      double aL =  Sr;
      double aR = -Sl;

      prim2cons( primL , Ul , r , 1.0 );
      prim2cons( primR , Ur , r , 1.0 );
      flux( primL , Fl , r , n );
      flux( primR , Fr , r , n );

      for( q=0 ; q<NUM_Q ; ++q ){

         Fstar = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
         Ustar = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );

         Flux[q] = Fstar - w*Ustar;

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
   double dAdt      = F->dA*dt;
   double dxL       = F->dxL;
   double dxR       = F->dxR;
   double r         = F->cm[0];

   double primL[NUM_Q];
   double primR[NUM_Q];

   double pL = cL->riph - .5*cL->dr;
   double pR = cR->riph - .5*cR->dr;
   double drL = r - pL;
   double drR = pR - r;

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + cL->grad[q]*dxL + cL->gradr[q]*drL;
      primR[q] = cR->prim[q] - cR->grad[q]*dxR - cR->gradr[q]*drR;
   }
   
   double Sl,Sr,Ss;
   double n[3] = {0.0,0.0,0.0};
   n[dim] = 1.0;

   vel( primL , primR , &Sl , &Sr , &Ss , n , r );

   double Fl[NUM_Q];
   double Fr[NUM_Q];
   double Ul[NUM_Q];
   double Ur[NUM_Q];
   
   double Flux[NUM_Q];

   if( 0. < Sl ){
      flux( primL , Fl , r , n ); 

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fl[q];
      }    
   }else if( 0. > Sr ){
      flux( primR , Fr , r , n ); 
     
      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fr[q];
      }
   }else{

      double aL =  Sr;
      double aR = -Sl;

      prim2cons( primL , Ul , r , 1.0 );
      prim2cons( primR , Ur , r , 1.0 );
      flux( primL , Fl , r , n );
      flux( primR , Fr , r , n );

      for( q=0 ; q<NUM_Q ; ++q ){

         Flux[q] = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );

      }

   }

   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dAdt;
      cR->cons[q] += Flux[q]*dAdt;
   }

}

