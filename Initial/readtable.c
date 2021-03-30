
#include "../paul.h"

static int NL = 0;
static double * rr  = NULL;
static double * rho = NULL;
static double * Pp  = NULL;
static double * vr  = NULL;
static double * Om  = NULL;
static double * Om_He  = NULL;
static double * Om_Ni  = NULL;

int countlines(char * filename){
   FILE *pFile = fopen(filename, "r");
   int lines=0;
   char c;
   while ((c = fgetc(pFile)) != EOF){
      if (c == '\n') ++lines;
   }
   fclose(pFile);
   return(lines);
}

int getTable( void ){
   int nL = countlines("Initial/initial.dat");
   rr  = (double *) malloc( nL*sizeof(double) );
   rho = (double *) malloc( nL*sizeof(double) );
   Pp  = (double *) malloc( nL*sizeof(double) );
   vr  = (double *) malloc( nL*sizeof(double) );
   Om  = (double *) malloc( nL*sizeof(double) );
   Om_He  = (double *) malloc( nL*sizeof(double) );
   Om_Ni  = (double *) malloc( nL*sizeof(double) );
   FILE * pFile = fopen("Initial/initial.dat","r");
   int l;
   for( l=0 ; l<nL ; ++l ){
      fscanf(pFile,"%lf %lf %lf %lf %lf %lf %lf\n",&(rr[l]),&(rho[l]),&(Pp[l]),&(vr[l]),&(Om[l]),&(Om_He[l]),&(Om_Ni[l]));
   }
   fclose(pFile);
   return(nL);
}

static double R_MIN;
static double EXP_ENERGY_B;
static double GAMMA_LAW;

void setICparams( struct domain * theDomain ){
   NL = getTable();
   R_MIN = theDomain->theParList.rmin;
   EXP_ENERGY_B = theDomain->theParList.Explosion_Energy;
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){

//   double theta_0 = 0.1;

   double r=x[0];
//r += 1e8;
   int l=0;
   while( rr[l] < r && l < NL-2 ) ++l;
   if( l==0 ) ++l;

   double rp = rr[l];
   double rm = rr[l-1];
   double drm = fabs(r-rm);
   double drp = fabs(rp-r);

   double P     = (Pp[l-1]*drp  + Pp[l]*drm )/(drp+drm);
   double rh    = (rho[l-1]*drp + rho[l]*drm)/(drp+drm);
   double V     = (vr[l-1]*drp  + vr[l]*drm )/(drp+drm);
   double X     = (Om[l-1]*drp  + Om[l]*drm )/(drp+drm);
   double X_He  = (Om_He[l-1]*drp  + Om_He[l]*drm )/(drp+drm);
   double X_Ni  = (Om_Ni[l-1]*drp  + Om_Ni[l]*drm )/(drp+drm);

   if( V < 0.0 ) V = 0.0;
   if( l==1 ) V = 0.0;
/*
   if( r < rr[1] ){
      rh = rho[1]*pow(r/rr[1],-6.);
      P  = Pp[1]*pow(r/rr[1],-6.);
//      V  = vr[1]*pow(r/rr[1],5.);
   }
*/
   if( l==NL-2 ){
      V = 0.0;
      rh = 1e18/r/r;
      P  = rh*Pp[NL-2]/rho[NL-2];
      //X = 1.0;
      X = Om[NL-2];
      X_He = Om_He[NL-2];
      X_Ni = Om_Ni[NL-2];
   }

//   if( X > 0.5 ) X = 0.0; else X = 1.0;

   //double P = rh*1e-8;
/* 
   double th = x[1];
   double w = exp(-pow(th/theta_0,30.));

   double rho0 = 1.0;
   double P0 = 1e-5;
   rh = rh*w + rho0*(1.-w);
   P  = P*w  + P0*(1.-w);
   V  = V*w;
*/
//   if( rh > 1e4 ) rh = 1e4;
//   if( P/rh > 1e17 ) P = 1e17*rh;
   double delta = (double)rand()/(double)RAND_MAX - .5;
   delta *= 1e-3;

   prim[RHO] = rh;
   prim[PPP] = P;
   prim[UU1] = V*(1.+delta);
   prim[UU2] = 0.0;
   if( NUM_C > 4 ) prim[UU3] = 0.0;

   if( NUM_N > 0 ){
//      double rtrs = 1e-3;
//      prim[NUM_C] = 0.0;
//      if( r > rtrs ) prim[NUM_C] = 1.0;
      prim[NUM_C] = X;
   }
   if( NUM_N > 1 ) prim[NUM_C+1] = X_He;
   if( NUM_N > 2 ) prim[NUM_C+2] = X_Ni;
}

void freeTable( void ){
   free(rr);
   free(rho);
   free(Pp);
   free(vr);
   free(Om);
}
