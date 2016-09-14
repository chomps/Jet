
#include <stdio.h>
#include "../paul.h"

double get_dL( double * a , double * b , double * c ){
   return(0.);
}

void prim2cons( double * , double * , double , double , double );
void cons2prim( double * , double * , double , double , double );

int main(void){

   double r  = 2.0;
   double dV = 3.0;

   double prim[NUM_Q];
   double cons[NUM_Q];
   double pnew[NUM_Q];
   
   prim[RHO] = 1.0;

   int q;
   srand(666);
   rand();
   double wiggle = .999;

   double wexp_min = -4.0;
   double wexp_max = 5.0;

   double pexp_min = -15.0;
   double pexp_max = 30.0;

   int i,j;
   int NW = 400;
   int NP = 400;
   for( i=0 ; i<NW ; ++i ){
      for( j=0 ; j<NP ; ++j ){

         double pexp = pexp_min + ((double)j+.5)*(pexp_max-pexp_min)/(double)NP;
         double wexp = wexp_min + ((double)i+.5)*(wexp_max-wexp_min)/(double)NW;
   
         double Pp = pow(10.,pexp);
         double W  = pow(10.,wexp);
         prim[PPP] = Pp;
         prim[UU1] = W;
         prim[UU2] = .01*W;

         for( q=0 ; q<NUM_Q ; ++q ){
            pnew[q] = prim[q];
            pnew[q] *= 1.0 + 2.*wiggle*( (double)rand()/(double)RAND_MAX - .5 );
         }

         prim2cons( prim , cons , r , dV );
         cons2prim( cons , pnew , r , dV );
   
         double err = 0.0;
         for( q=0 ; q<NUM_Q ; ++q ){
            err += pow( (prim[q] - pnew[q])/prim[q] , 2. );
         }
         err = sqrt(err);

         if( err>1e-6 ) printf("%e\t%e\t%e\n",W,Pp,err);
      }
   }

   return(0);
}

