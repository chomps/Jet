
#include "paul.h"

double cart_r( double x , double rmin , double rmax ){
   return( rmin + x*(rmax-rmin) );
}

double log_r( double x , double rmin , double rmax ){
   return( rmin*pow(rmax/rmin,x) );
}

double hybrid_r( double x , double rmin , double rmax , double R ){
   return( R*pow(rmax/R,x) + rmin-R + (R-rmin)*x );
}

double target_x( double x , double x0 , double weight ){
   x0 = 2.*x0-1.;
   x = 2.*x-1.;
   double y = (x-3.*x0)*(x*x-1.)/(3.*x0*x0+1.) + x;
   y = .5*(y+1.);
   x = .5*(x+1.);
   return( weight*y + (1.-weight)*x );
}

void start_clock( struct domain * );

void setupDomain( struct domain * theDomain ){

   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   theDomain->theCells = (struct cell **) malloc( Np*Nt*sizeof(struct cell *) );
   int j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         //theDomain->theCells[j+Nt*k] = (struct cell *) malloc( Nr[j+Nt*k]*sizeof(struct cell) );
         theDomain->theCells[j+Nt*k] = (struct cell *) calloc( Nr[j+Nt*k] , sizeof(struct cell) );
      }
   }
   int i;
   double R_MIN = theDomain->theParList.rmin;
   double R_MAX = theDomain->theParList.rmax;
   double LogWeight = theDomain->theParList.LogZoning;
   double X0 = theDomain->theParList.Target_X;
   double W0 = theDomain->theParList.Target_W;
   for( j=0 ; j<Nt ; ++j ){
      int offset=1;
      if( j%2 == 0 ) offset=-offset;
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         for( i=0 ; i<Nr[jk] ; ++i ){
            double x = (double)(i+1)/(double)Nr[jk];
            if( i!=0 && i!=Nr[jk]-1 ) x += 0.5*.125*(double)offset/(double)Nr[jk];
            x = target_x( x , X0 , W0 );
            double rx = cart_r(x,R_MIN,R_MAX);//R_MIN + x*(R_MAX-R_MIN);
            double rl = log_r(x,R_MIN,R_MAX);//R_MAX*pow(R_MAX/R_MIN,x-1.);
            double rp = rx;
            if( LogWeight > 0.5 ) rp = rl;
            //rp = hybrid_r( x , R_MIN , R_MAX , 0.1 );
            theDomain->theCells[j+Nt*k][i].riph = rp;//LogWeight*rl + (1.-LogWeight)*rx;
            theDomain->theCells[j+Nt*k][i].wiph = 0.0;
         }
      }
   }

   theDomain->g_point_mass = theDomain->theParList.PointMass*2.0e33;

   theDomain->t      = theDomain->theParList.t_min;
   theDomain->t_init = theDomain->theParList.t_min;
   theDomain->t_fin  = theDomain->theParList.t_max;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;

   theDomain->final_step = 0;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;

   theDomain->count_steps = 0;
   start_clock( theDomain );

}

void clear_cell( struct cell * c ){
   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      c->prim[q]   = 0.0;
      c->cons[q]   = 0.0;
      c->RKcons[q] = 0.0;
      c->grad[q]   = 0.0;
      c->gradr[q]  = 0.0;
   }
   c->riph = 0.0;
   c->RKriph = 0.0;
   c->dr = 0.0;
   c->wiph = 0.0;
}

void freeDomain( struct domain * theDomain ){
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         free( theDomain->theCells[j+Nt*k] );
      }
   }
   free( theDomain->theCells );
   free( theDomain->Nr );
   theDomain->t_jph--;
   free( theDomain->t_jph );
   theDomain->p_kph--;
   free( theDomain->p_kph );
}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final=0;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }

   if( theDomain->rank==0 ){
      FILE * abort = NULL;
      abort = fopen("abort","r");
      if( abort ){ final = 1; fclose(abort); }
   }

   MPI_Allreduce( MPI_IN_PLACE , &final , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   if( final ) theDomain->final_step = 1;

}

void report( struct domain * , double );
void snapshot( struct domain * , char * );
void output( struct domain * , char * );

void possiblyOutput( struct domain * theDomain , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   double Nsnp = theDomain->N_snp;
   double Nchk = theDomain->N_chk;
   int LogOut = theDomain->theParList.Out_LogTime;
   int n0;

   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      //longandshort( &theDomain , &L , &S , &iL , &iS , theDomain.theCells[0] , 0 , 0 );
      report( theDomain , t );
      int Nr = theDomain->Nr[0];
      double rout = theDomain->theCells[0][Nr-1].riph;
      if( theDomain->rank==0 ) printf("t = %.4e Nr = %d Rout = %.2e\n",t,Nr,rout);
   }
   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         if(theDomain->rank==0) printf("Creating Checkpoint #%04d...\n",n0);
         sprintf(filename,"checkpoint_%04d",n0);
         output( theDomain , filename );
      }else{
         if(theDomain->rank==0) printf("Creating Final Checkpoint...\n");
         output( theDomain , "output" );
      }
   }
   n0 = (int)( t*Nsnp/t_fin );
   if( LogOut ) n0 = (int)( Nsnp*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nsnp < n0 && Nsnp>0) || override ){
      theDomain->nsnp = n0;
      char filename[256];
      if(!override) sprintf( filename , "snapshot_%04d" , n0 );
      else sprintf( filename , "snapshot" );
      snapshot( theDomain , filename );
   }
//printf("rank = %d t = %e n0 = %d\n",theDomain->rank,t,n0);

}
