#include "../paul.h"

static double GAMMA_LAW = 0.0;
static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;

void setHydroParams( struct domain * theDomain ){
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
}

double get_vr( double * prim ){
   return( prim[UU1] );
}

double get_entropy( double * prim ){
   return( log( prim[PPP] / pow( prim[RHO] , GAMMA_LAW ) ) );
}

double get_pot( double r );

void prim2cons( double * prim , double * cons , double r , double th , double dV ){
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[UU1];
   double vt  = prim[UU2];
   double v2 = vr*vr + vt*vt;
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   double phi = get_pot( r );

   cons[DEN] = rho*dV;
   cons[SS1] = rho*vr*dV;
   cons[SS2] = r*rho*vt*dV;
   cons[TAU] = (.5*rho*v2 + rhoe + rho*phi )*dV;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
      cons[q] = cons[DEN]*prim[q];
   }
}

void cons2prim( double * cons , double * prim , double r , double th , double dV ){

   double rho = cons[DEN]/dV;
   double Sr  = cons[SS1]/dV;
   double St  = cons[SS2]/dV/r;
   double E   = cons[TAU]/dV;

   double vr = Sr/rho;
   double vt = St/rho;
   double v2 = vr*vr + vt*vt;

   double phi = get_pot( r );

   double rhoe = E - .5*rho*v2 - rho*phi;

   double gam = GAMMA_LAW;
   double Pp = (gam-1.)*rhoe;

   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   if( Pp < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;


   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = vr;
   prim[UU2] = vt;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
      prim[q] = cons[q]/cons[DEN];
   }
}

void getUstar( double * prim , double * Ustar , double r , double th , double Sk , double Ss , double * n ){

   double rho = prim[RHO];
   double vr  = prim[UU1];
   double vt  = prim[UU2];
   double Pp  = prim[PPP];
   double v2  = vr*vr+vt*vt;

   double vn = vr*n[0] + vt*n[1];
   double gam = GAMMA_LAW;

   double phi = get_pot( r );
   double rhoe = Pp/(gam-1.);

   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   Ustar[DEN] = rhostar;
   Ustar[SS1] = rhostar*( vr + (Ss-vn)*n[0] );
   Ustar[SS2] = r*rhostar*( vt + (Ss-vn)*n[1] );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vn) + rhostar*phi + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DEN];
   }

}

void flux( double * prim , double * flux , double r , double th , double * n ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[UU1];
   double vt  = prim[UU2];
   double v2  = vr*vr + vt*vt;
   double vn  = vr*n[0] + vt*n[1];
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   double phi = get_pot( r );
 
   flux[DEN] = rho*vn;
   flux[SS1] = rho*vr*vn + Pp*n[0];
   flux[SS2] = r*( rho*vt*vn + Pp*n[1] );
   flux[TAU] = (.5*rho*v2 + rhoe + rho*phi + Pp)*vn;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
      flux[q] = flux[DEN]*prim[q];
   }
}

void source( double * prim , double * cons , double * xp , double * xm , double dVdt ){
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double rp = xp[0];
   double rm = xm[0];
   double r  = .5*(rp+rm);
   double th = .5*(xp[1]+xm[1]);
   double vt  = prim[UU2];
   double r2 = (rp*rp+rm*rm+rp*rm)/3.;
   cons[SS1] += (2.*Pp + rho*vt*vt)*(r/r2)*dVdt;
   cons[SS2] += Pp*(cos(th)/sin(th))*dVdt;//*(r*r/r2);
}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss , double * n , double r , double th ){
   
   double gam = GAMMA_LAW;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vr1  = prim1[UU1];
   double vt1  = prim1[UU2];
   double vn1  = vr1*n[0]+vt1*n[1];

   double cs1 = sqrt(fabs(gam*P1/rho1));

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vr2  = prim2[UU1];
   double vt2  = prim2[UU2];
   double vn2  = vr2*n[0]+vt2*n[1];

   double cs2 = sqrt(fabs(gam*P2/rho2));

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;
   
}

double get_dL( double * , double * , int );

double mindt( double * prim , double w , double * xp , double * xm ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[UU1];
   double vt  = prim[UU2];//*.5*(xp[0]+xm[0]);
   double gam = GAMMA_LAW;

   double cs = sqrt(fabs(gam*Pp/rho));

   double maxvr = cs + fabs( vr - w );
   double maxvt = cs + fabs( vt );
   double dtr = get_dL(xp,xm,0)/maxvr;
   double dtt = get_dL(xp,xm,1)/maxvt;
   double dt = dtr;
   if( dt > dtt ) dt = dtt;

   return( dt );

}
