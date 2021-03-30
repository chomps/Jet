
#include "paul.h"

int mpiSetup( struct domain * , int , char *[] );
void gridSetup( struct domain * );
void timestep( struct domain * , double );
void setupcells( struct domain * );
void regrid( struct domain * );
void exchangeData( struct domain * , int );
double getmindt( struct domain * );

void setICparams( struct domain * );
void setHydroParams( struct domain * );
void setNozzleParams( struct domain * );

void read_par_file( struct domain * );

void setupDomain( struct domain * );
void freeDomain( struct domain * );
void check_dt( struct domain * , double * );
void possiblyOutput( struct domain * , int );
void generate_log( struct domain * );

int main( int argc , char * argv[] ){
   
   MPI_Init(&argc,&argv);
   struct domain theDomain = {0};
   read_par_file( &theDomain );
   
   int error = mpiSetup(&theDomain,argc,argv);
   if( error==1 ) return(0);

   if(theDomain.rank==0) remove("abort");

   gridSetup( &theDomain );   

   setupDomain( &theDomain );
   setICparams( &theDomain );
   setHydroParams( &theDomain );
   setNozzleParams( &theDomain );
   setupcells( &theDomain );

   if( theDomain.theParList.Initial_Regrid && !(theDomain.theParList.restart_flag) ) regrid( &theDomain );

   if( theDomain.Nt > 1 ) exchangeData( &theDomain , 0 );
   if( theDomain.Np > 1 ) exchangeData( &theDomain , 1 );

   if( theDomain.rank==0 && !(theDomain.theParList.restart_flag) ){
      FILE * rFile = fopen("report.dat","w");
      fclose(rFile);
   } 


   while( !(theDomain.final_step) ){

      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      possiblyOutput( &theDomain , 0 );
      timestep( &theDomain , dt );

   }

   possiblyOutput( &theDomain , 1 );
   generate_log( &theDomain );
   MPI_Barrier(theDomain.theComm);
   MPI_Comm_free(&(theDomain.theComm));
   MPI_Finalize();
   freeDomain( &theDomain );
   return(0);

}

