
#include "paul.h"

double get_dA( double * , double * , int );

int get_num_tpFaces( int Nt , int Np , int dim ){
   if( dim==1 ) return( (Nt-1)*Np );
   else return( Nt*Np );
}

void addFace( struct face * theFaces , int n , struct cell * cL , struct cell * cR , double dxL , double dxR , double * xp , double * xm , int dim ){
   int d;
   for( d=0 ; d<3 ; ++d ) theFaces[n].cm[d] = .5*(xp[d]+xm[d]);
   double rp = xp[0];
   double rm = xm[0];
   double r2 = (rp*rp+rm*rm+rp*rm)/3.;
   theFaces[n].cm[0] = r2/(.5*(rp+rm));
   theFaces[n].L   = cL;
   theFaces[n].R   = cR;
   theFaces[n].dxL = dxL;
   theFaces[n].dxR = dxR;
   theFaces[n].dr  = xp[0]-xm[0];
   theFaces[n].dA  = get_dA(xp,xm,dim);
}

void buildfaces( struct domain * theDomain , struct face * theFaces , int * ntj , int dim , int mode ){
  
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int i,j,k; 
   int n=0;
 
   int Ntmax = Nt-1;
   if( dim==2 ) Ntmax = Nt;

   for( j=0 ; j<Ntmax ; ++j ){
      int jp = j+1;
      for( k=0 ; k<Np ; ++k ){
         int JK = j+Ntmax*k;
         if( mode == 0 ) ntj[JK] = n;

         int kp = k+1;
         if( kp == Np ) kp = 0;
         int jk  = j  + Nt*k;
         int jkp = jp + Nt*k;
         if( dim==2 ) jkp = j + Nt*kp;

         double dxL = .5*(t_jph[j]  - t_jph[j-1]);
         double dxR = .5*(t_jph[jp] - t_jph[j]  );
         if( dim==2 ){
            dxL = .5*(p_kph[k]  - p_kph[k-1]);
            dxR = .5*(p_kph[kp] - p_kph[k]  );
            if( dxR < 0. ) dxR += M_PI;
         }
         double xp[3] = {0.0,t_jph[j],p_kph[k  ]};
         double xm[3] = {0.0,t_jph[j],p_kph[k-1]};
         if( dim==2 ){ xm[1] = t_jph[j-1] ; xm[2] = p_kph[k] ; }

         int ip=0;
         for( i=0 ; i<Nr[jk] ; ++i ){
            struct cell * cL = &(theCells[jk ][i] );
            struct cell * cR = &(theCells[jkp][ip]);
            //First figure out if cell+ covers all of cell-, if so create one face out of cell-, and move on to next i.
            if( cR->riph > cL->riph ){
               xp[0] = cL->riph;
               xm[0] = cL->riph-cL->dr;
               if( mode==1 ) addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim );
               ++n;
            }else{
            //Otherwise, three steps:
               //Step A: face formed out of beginning of cell- and end of cell+. ++ip;
               xp[0] = cR->riph;
               xm[0] = cL->riph-cL->dr;
               if( mode==1 ) addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim );
               ++n;

               ++ip;
               if( ip < Nr[jkp] ) cR = &(theCells[jkp][ip]);
               double dr = cL->riph - cR->riph;
               while( dr > 0.0 && ip < Nr[jkp] ){
                  //Step B: (optional) all faces formed out of part of cell- and all of cell+. ++ip;
                  xp[0] = cR->riph;
                  xm[0] = cR->riph-cR->dr;
                  if( mode==1 ) addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim );
                  ++n;
               
                  ++ip;
                  if( ip < Nr[jkp] ) cR = &(theCells[jkp][ip]);
                  dr = cL->riph - cR->riph;
               }
               if( ip < Nr[jkp] ){
                  //Step C: face formed out of end of cell- and beginning of cell+.
                  xp[0] = cL->riph;
                  xm[0] = cR->riph-cR->dr;
                  if( mode==1 ) addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim );
                  ++n;
               }
            }
         }
      }
   }
   if( mode==0 ){
      int NN = get_num_tpFaces( Nt , Np , dim );
      ntj[NN] = n;
   }
}


