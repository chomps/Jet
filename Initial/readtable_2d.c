
#include "../paul.h"

static int Nx = 0;
static double * rho = NULL;
static double xmax;
static double Mass;

void getTable( void ){
   FILE * pFile = fopen("Initial/table.dat","r");
   fscanf(pFile,"%d %lf %lf\n",&Nx,&Mass,&xmax);
   rho = (double *) malloc( Nx*Nx*sizeof(double) );
   int l;
   for( l=0 ; l<Nx*Nx ; ++l ){
      int i,j;
      double rho0;
      fscanf(pFile,"%d %d %lf\n",&i,&j,&rho0);
      rho[i*Nx+j] = rho0;
   }
   fclose(pFile);
}

void initial( double * prim , double * x ){

   double Ein = 28.0;       //x 10^51 erg
   double Mstar = 4.77;//8.0;    //Mass of star in solar masses
   double Rexp = R_MIN;

   double Msolc2 = 1790.0; //Solar mass * c^2 / 10^51 erg

   double EoverM = Ein/Mstar/Msolc2;

   double rho0 = 1.0;
   double rhoc = 0.003*rho0;
   double rho_min = 1e-4*rho0;
   double p_min = 1e-8;

   double r  = x[0];
   double th = x[1];
   if( Nx==0 ) getTable();

   double E = EoverM*Mass;
   double e0 = E/pow(sqrt(M_PI)*Rexp,3.);
   double ex = e0*exp(-r*r/Rexp/Rexp);
   double Pexp = (GAMMA_LAW-1.)*ex;

   double s = r*sin(th);
   double z = r*cos(th);

   double ix,jy;
   ix = (double)Nx/xmax*s;
   jy = (double)Nx/xmax*fabs(z);

   int i0 = (int)ix;
   if( i0 > Nx-2 ) i0 = Nx-2;
   if( i0 < 0    ) i0 = 0;
   int i1 = i0+1;
   int j0 = (int)jy;
   if( j0 > Nx-2 ) j0 = Nx-2;
   if( j0 < 0    ) j0 = 0;
   int j1 = j0+1;

   double xx = ix-i0;
   if( xx > 1. ) xx = 1.;
   if( xx < 0. ) xx = 0.;
   double yy = jy-j0;
   if( yy > 1. ) yy = 1.;
   if( yy < 0. ) yy = 0.;

   double fa = rho[i0*Nx+j0];
   double fb = rho[i0*Nx+j1];
   double fc = rho[i1*Nx+j1];
   double fd = rho[i1*Nx+j0];

   double rho_interp = xx*yy*( (fa+fc) - (fb+fd) ) + xx*(fd-fa) + yy*(fb-fa) + fa;

   double rho = rho0*rho_interp;

   if( rho < rhoc ){
      double XX = rho/rhoc-.5;
      if( XX < 0.0 ) XX = 0.0;
      rho = rhoc * sqrt(2.*XX);
   }

   //if( rho < rho_min ) rho = rho_min;

   double R = 1.0;
   double Msol = Mass/Mstar;
   double Msol_py = 1.49e-8*Msol/R;
   double Mdot = 1e-5*Msol_py;
   double vwind = .012;
   double wind = Mdot/4./M_PI/vwind/r/r;

   rho += wind;

   double Pp = Pexp + rho*p_min;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = 0.0;
   prim[UU2] = 0.0;

}

void freeTable( void ){
   free(rho);
}
