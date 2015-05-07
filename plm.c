
#include "paul.h"

double minmod( double a , double b , double c ){
   double m = a;
   if( a*b < 0.0 ) m = 0.0;
   if( fabs(b) < fabs(m) ) m = b;
   if( b*c < 0.0 ) m = 0.0;
   if( fabs(c) < fabs(m) ) m = c;
   return(m);
}

void plm_r( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double PLM = theDomain->theParList.PLM;
   int i,j,k,q;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         for( i=0 ; i<Nr[jk] ; ++i ){
            struct cell * c  = &(theCells[jk][i]);
            if( i==0 || i==Nr[jk]-1 ){
               for( q=0 ; q<NUM_Q ; ++q ){
                  c->gradr[q] = 0.0;
               }
            }else{
               struct cell * cL = &(theCells[jk][i-1]);
               struct cell * cR = &(theCells[jk][i+1]);
               double drL = cL->dr;
               double drC = c->dr;
               double drR = cR->dr;
               for( q=0 ; q<NUM_Q ; ++q ){
                  double pL = cL->prim[q];
                  double pC = c->prim[q];
                  double pR = cR->prim[q];
                  double sL = pC - pL;
                  sL /= .5*( drC + drL );
                  double sR = pR - pC;
                  sR /= .5*( drR + drC );
                  double sC = pR - pL;
                  sC /= .5*( drL + drR ) + drC;
                  sL *= PLM;
                  sR *= PLM;
                  c->gradr[q] = minmod( sL , sC , sR );
               }
            }
         }
      }
   }
}

double get_dA( double * , double * , int );

void plm_trans( struct domain * theDomain , struct face * theFaces , int Nf , int dim ){

   struct cell ** theCells = theDomain->theCells;
   int Np = theDomain->Np;
   int Nt = theDomain->Nt;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double PLM = theDomain->theParList.PLM;
   int i,j,k,q;

   //Clear gradients
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         for( i=0 ; i<Nr[jk] ; ++i ){
            for( q=0 ; q<NUM_Q ; ++q ){
               theCells[jk][i].grad[q] = 0.0;
            }
         }
      }
   }

   //Add weighted slopes
   int n;
   for( n=0 ; n<Nf ; ++n ){
      struct face * f  = &( theFaces[n] );
      double r = f->cm[0];
      struct cell * cL = f->L;
      struct cell * cR = f->R;
      double dxL = f->dxL;
      double dxR = f->dxR;
      double rL = cL->riph - .5*cL->dr;
      double rR = cR->riph - .5*cR->dr;
      double drL = r - rL;
      double drR = rR - r;
      double dA  = f->dA;
      for( q=0 ; q<NUM_Q ; ++q ){
         double WL = cL->prim[q] + drL*cL->gradr[q];
         double WR = cR->prim[q] - drR*cR->gradr[q];

         double S = (WR-WL)/(dxR+dxL);
         cL->grad[q] += S*dA;
         cR->grad[q] += S*dA;
      }
   }

   //Divide by total weight
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         int j_outside = 0;
         if( j==0 || j==Nt-1 ) j_outside = 1;
         for( i=0 ; i<Nr[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double rp = c->riph;
            double rm = c->riph-c->dr;
            double xp[3] = {rp,t_jph[j  ],p_kph[k  ]};
            double xm[3] = {rm,t_jph[j-1],p_kph[k-1]};
            if( dim==1 ) xm[1] = t_jph[j];
            else xm[2] = p_kph[k];
            double dA1 = get_dA(xp,xm,dim);
            if( dim==1 ){
               xp[1] = t_jph[j-1];
               xm[1] = t_jph[j-1];
            }else{
               xp[2] = p_kph[k-1];
               xm[2] = p_kph[k-1];
            }
            double dA2 = get_dA(xp,xm,dim);
            double dAtot = dA1 + dA2;
            for( q=0 ; q<NUM_Q ; ++q ){
               c->grad[q] /= dAtot;
               if( j_outside ) c->grad[q] = 0.0;
            }    
         }    
      }    
   }

   //Slope Limiting
   for( n=0 ; n<Nf ; ++n ){
      struct face * f  = &( theFaces[n] );
      double r = f->cm[0];
      struct cell * cL = f->L;
      struct cell * cR = f->R;
      double dxL = f->dxL;
      double dxR = f->dxR;
      double rL = cL->riph - .5*cL->dr;
      double rR = cR->riph - .5*cR->dr;
      double drL = r - rL;
      double drR = rR - r;
      for( q=0 ; q<NUM_Q ; ++q ){
         double WL = cL->prim[q] + drL*cL->gradr[q];
         double WR = cR->prim[q] - drR*cR->gradr[q];

         double S = (WR-WL)/(dxR+dxL);
         double SL = cL->grad[q];
         double SR = cR->grad[q];
         if( S*SL < 0.0 ){
            cL->grad[q] = 0.0; 
         }else if( fabs(PLM*S) < fabs(SL) ){
            cL->grad[q] = PLM*S;
         }
         if( S*SR < 0.0 ){
            cR->grad[q] = 0.0; 
         }else if( fabs(PLM*S) < fabs(SR) ){
            cR->grad[q] = PLM*S;
         }    
      }    
   }
}


