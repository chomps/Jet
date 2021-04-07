
#include "paul.h"

//#define G_CONST     6.673e-8
#define G_CONST     0.0
//#define G_NUM_R     1000 //This should really be NUM_R
//#define G_EPS       0.0

static int G_NUM_R=0;
static double * Menc = NULL;
static double * Pot  = NULL;
//static double Menc[G_NUM_R];

void clear_menc( double g_point_mass ){
   int i;
   for( i=0 ; i<G_NUM_R ; ++i ){
      Menc[i] = 0.0;
   }
   Menc[0] = g_point_mass;
}

double xtor( double x , double rmin , double rmax ){
   return( rmin*pow(rmax/rmin,x) );
}

double rtox( double r , double rmin , double rmax ){
   return( log( r/rmin ) / log( rmax/rmin ) );
}

double menc_force( double x , double r ){

   int ir = (int)(x*G_NUM_R);
   if( ir<0 ) ir=0;
   if( ir>G_NUM_R-1 ) ir = G_NUM_R-1;
   double dx = x*((double)G_NUM_R)-(double)ir;
   if( dx<0.0 ) dx=0.0;
   if( dx>1.0 ) dx=1.0;
   double M0 = Menc[ir];
   double M1;
   if(ir<G_NUM_R-1) M1 = Menc[ir+1]; else M1 = Menc[ir];
   double M = M0*(1.-dx) + M1*dx;

//
double rhoc = 1./4./M_PI/sqrt(3.);
   M = 4.*M_PI/3.*rhoc*r*r*r/pow(1. + r*r/3.,1.5);
//

   double f = -G_CONST*M/r/r;

   return(f);

}

double get_pot( double r ){
   double rhoc = 1./4./M_PI/sqrt(3.);
   double phi = -4.*M_PI*G_CONST*rhoc/sqrt( 1. + r*r/3. );
phi = 0.0;

   return( phi );
}

void grav_src( double * cons , double dt , double r , double x ){

   double m  = cons[DEN];
   double Sr = cons[SS1];
   
   double f = menc_force( x , r );
/* 
   double GM = 6.67e-8*5.45*2e33;
   double R = 2.47*7e10;
   double f = -GM/r/r;
   if( r<R ) f = 0.0;
*/
   //printf("F=%e\n",m*f*dt);

   cons[SS1] += m*f*dt;
//   cons[TAU] += Sr*f*dt;

}

void acct_menc( struct cell * thisCell , double rmin , double rmax ){

   double r = thisCell->riph - .5*thisCell->dr;
   double x = rtox( r , rmin , rmax );
   int ir = (int) (x*G_NUM_R);
   if( ir<0 ) ir=0;

   if( ir<G_NUM_R-1 ) Menc[ir] += thisCell->cons[DEN];

}

void aggregate_mass( void ){
   int i;
   double Mtot[G_NUM_R];
   double M=0.0;
   for( i=0 ; i<G_NUM_R ; ++i ){
      M += Menc[i];
      Mtot[i] = M;
   }

   MPI_Allreduce( Mtot , Menc , G_NUM_R , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
}

void calculate_mass( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;

   clear_menc( theDomain->g_point_mass );

   int j,k,i;

   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         double rmin = theCells[jk][0].riph;
         double rmax = theCells[jk][Nr[jk]-1].riph;
         for( i=1 ; i<Nr[jk] ; ++i ){
            acct_menc( &(theCells[jk][i]) , rmin , rmax );
         }
      }
   }

   aggregate_mass();

}

void gravity_setup( struct domain * theDomain ){

   if( !G_NUM_R ) G_NUM_R = theDomain->theParList.Num_R;
   if( !Menc ) Menc = (double *) malloc(G_NUM_R*sizeof(double));

   calculate_mass( theDomain );

}

void gravity_outputMass( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr[0];

   double Rmin = theCells[0][ 0  ].riph;
   double Rmax = theCells[0][Nr-1].riph;
   FILE * pFile = fopen("mass.dat","w");
   int i;
   for( i=0 ; i<G_NUM_R ; ++i ){
      double x = (double)i/(double)G_NUM_R;
      double r = xtor(x,Rmin,Rmax);
      fprintf(pFile,"%e %e\n",r,Menc[i]);
   }
   fclose(pFile);

}

void gravity_addsrc( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;

   int i,jk;
   for( jk=0 ; jk<Nt*Np ; ++jk ){
      double rmin = theCells[jk][0].riph;
      double rmax = theCells[jk][Nr[jk]-1].riph;
      for( i=0 ; i<Nr[jk]-1 ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         double r = c->riph-.5*c->dr;
         double x = rtox( r , rmin , rmax );
         grav_src( c->cons , dt , r , x );
      }
   }

   free(Menc);
   Menc = NULL;
}

