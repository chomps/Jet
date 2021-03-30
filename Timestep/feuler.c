#include "../paul.h"

void onestep( struct domain * , double , double , int , double );

void timestep( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;

   int i,jk;
   for( jk=0 ; jk<Nt*Np ; ++jk ){
      for( i=0 ; i<Nr[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
         c->RKriph = c->riph; 
      }    
   }
   onestep( theDomain , 0.0 , dt , 1 , dt );

   theDomain->t += dt;
   theDomain->count_steps++;

}
