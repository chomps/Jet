#include "../paul.h"

void prim2cons( double * prim , double * cons , double r , double dV ){
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double v2 = vx*vx + vy*vy;
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);
   cons[DEN] = rho*dV;
   cons[SS1] = rho*vx*dV;
   cons[SS2] = rho*vy*dV;
   cons[TAU] = (.5*rho*v2 + rhoe)*dV;
}

void cons2prim( double * cons , double * prim , double r , double dV ){
   double rho = cons[DEN]/dV;
   double Sx  = cons[SS1]/dV;
   double Sy  = cons[SS2]/dV;
   double E   = cons[TAU]/dV;

   double vx = Sx/rho;
   double vy = Sy/rho;
   double v2 = vx*vx + vy*vy;
   double rhoe = E - .5*rho*v2;
   double gam = GAMMA_LAW;
   double Pp = (gam-1.)*rhoe;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = vx;
   prim[UU2] = vy;

}

void getUstar( double * prim , double * Ustar , double r , double Sk , double Ss , double * n ){

   double rho = prim[RHO];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double Pp  = prim[PPP];
   double v2 = vx*vx+vy*vy;

   double vn = vx*n[0] + vy*n[1];
   double gam = GAMMA_LAW;

   double rhoe = Pp/(gam-1.);

   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   Ustar[DEN] = rhostar;
   Ustar[SS1] = rhostar*( vx + (Ss-vn)*n[0] );
   Ustar[SS2] = rhostar*( vy + (Ss-vn)*n[1] );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vn) + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DEN];
   }

}

void flux( double * prim , double * flux , double r , double * n ){
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double v2  = vx*vx + vy*vy;
   double vn  = vx*n[0] + vy*n[1];
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);
   
   flux[DEN] = rho*vn;
   flux[SS1] = rho*vx*vn + Pp*n[0];
   flux[SS2] = rho*vy*vn + Pp*n[1];
   flux[TAU] = (.5*rho*v2 + rhoe + Pp)*vn;
}

void source( double * prim , double * cons , double * x , double dVdt ){
   //Silence is golden.
}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss , double * n , double r ){
   
   double gam = GAMMA_LAW;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vx1  = prim1[UU1];
   double vy1  = prim1[UU2];
   double vn1  = vx1*n[0]+vy1*n[1];

   double cs1 = sqrt(gam*P1/rho1);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vx2  = prim2[UU1];
   double vy2  = prim2[UU2];
   double vn2  = vx2*n[0]+vy2*n[1];

   double cs2 = sqrt(gam*P2/rho2);

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
   double vx  = prim[UU1];
   double vy  = prim[UU2];
   double gam = GAMMA_LAW;

   double cs = sqrt( gam*Pp/rho );

   double maxvx = cs + fabs( vx - w );
   double maxvy = cs + fabs( vy );
   double dtx = get_dL(xp,xm,0)/maxvx;
   double dty = get_dL(xp,xm,1)/maxvy;
   double dt = dtx;
   if( dt > dty ) dt = dty;

   return( dt );

}
