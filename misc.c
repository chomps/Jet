#include "paul.h"
#include <string.h>

double get_dA( double * , double * , int );
double get_dV( double * , double * );

double mindt( double * , double , double * , double * );

double getmindt( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;

   double dt = 1e100;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         for( i=1 ; i<Nr[j+Nt*k]-1 ; ++i ){
            struct cell * c = &(theCells[j+Nt*k][i]);
            double xp[3] = {c->riph      ,t_jph[j  ],p_kph[k  ]};
            double xm[3] = {c->riph-c->dr,t_jph[j-1],p_kph[k-1]};
            double wm = 0.0; 
            if( i>0 ) wm = theCells[j+Nt*k][i-1].wiph;
            double wp = c->wiph;
            double w = .5*(wm+wp);
            double dt_temp = mindt( c->prim , w , xp , xm );
            if( dt > dt_temp ) dt = dt_temp;
         }
      }
   }
   dt *= theDomain->theParList.CFL; 
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , theDomain->theComm );
   return( dt );
}

void initial( double * , double * );
void prim2cons( double * , double * , double , double , double );
void cons2prim( double * , double * , double , double , double );
void restart( struct domain * );
void exchangeData( struct domain * , int );

void setupcells( struct domain * theDomain ){

   int restart_flag = theDomain->theParList.restart_flag;
   if( restart_flag ){
      restart( theDomain );
if( theDomain->rank==0 ){ printf("Find the segfault 1...\n"); sleep(1); }
      if( theDomain->Np > 1 ) exchangeData( theDomain , 1 );
      if( theDomain->Nt > 1 ) exchangeData( theDomain , 0 );
if( theDomain->rank==0 ){ printf("Find the segfault 2...\n"); sleep(1); }
   }

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;

   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         for( i=0 ; i<Nr[j+Nt*k] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double rm = 0.0;
            if( i>0 ) rm = theCells[jk][i-1].riph;
            double rp = c->riph;
            c->dr = rp-rm;
            c->wiph = 0.0;
            double r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
            double xp[3] = {rp,t_jph[j]  ,p_kph[k]  };
            double xm[3] = {rm,t_jph[j-1],p_kph[k-1]};
            double dV = get_dV( xp , xm );
            double x[3] = {.5*(rp+rm),.5*(t_jph[j]+t_jph[j-1]),.5*(p_kph[k]+p_kph[k-1])};
            if( !restart_flag ){
               int icN = theDomain->theParList.Initial_Cons;
               if( !icN ){ 
                  initial( c->prim , x );
                  prim2cons( c->prim , c->cons , r , x[1] , dV );
               }else{
                  double consTot[NUM_Q];
                  int q;
                  for( q=0 ; q<NUM_Q ; ++q ) consTot[q] = 0.0;
                  int icn;
                  for( icn=0 ; icn<icN ; ++icn ){
                     double rmn = ( ((double)icn   )/(double)icN )*(rp-rm)+rm;
                     double rpn = ( ((double)icn+1.)/(double)icN )*(rp-rm)+rm;
                     double rn = (2./3.)*(rpn*rpn*rpn-rmn*rmn*rmn)/(rpn*rpn-rmn*rmn);
                     double xn[3] = {.5*(rpn+rmn),.5*(t_jph[j]+t_jph[j-1]),.5*(p_kph[k]+p_kph[k-1])};
                     double xpn[3] = {rpn,t_jph[j]  ,p_kph[k]  };
                     double xmn[3] = {rmn,t_jph[j-1],p_kph[k-1]};
                     double dVn = get_dV( xpn , xmn );
                     initial( c->prim , xn );
                     prim2cons( c->prim , c->cons , rn , xn[1] , dVn );
                     for( q=0 ; q<NUM_Q ; ++q ) consTot[q] += c->cons[q];
                  }
                  for( q=0 ; q<NUM_Q ; ++q ) c->cons[q] = consTot[q];
               }
            }else{
               prim2cons( c->prim , c->cons , r , x[1] , dV );
            }
            cons2prim( c->cons , c->prim , r , x[1] , dV );
         }
      }
   }
}

void clear_w( struct domain * theDomain ){
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         for( i=0 ; i<Nr[j+Nt*k] ; ++i ){
            theDomain->theCells[j+Nt*k][i].wiph = 0.0;
         }
      }
   }
}

double get_vr( double * );

void set_wcell( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;

   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         theCells[jk][0].wiph = 0.0;
         theCells[jk][Nr[jk]-1].wiph = 0.0;
         for( i=1 ; i<Nr[jk]-2 ; ++i ){
            struct cell * cL = &(theCells[jk][i  ]);  
            struct cell * cR = &(theCells[jk][i+1]);
            double wL = get_vr( cL->prim );
            double wR = get_vr( cR->prim );
            cL->wiph = .5*(wL+wR);
         }
      }    
   }
}

void initial( double * , double * );
void clear_cell( struct cell * );

void move_BCs( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double ShockPos = theDomain->theParList.ShockPos;

   int homologous_BCs = 0;
//   double r_outer = theCells[0][Nr[0]-1].riph + theCells[0][Nr[0]-1].dr;

   double max_vel = 0.0;
   int i,j,k; 
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         for( i=1 ; i<Nr[jk]-1 ; ++i ){
            double ww = get_vr( theCells[jk][i].prim);
            if( max_vel < ww ) max_vel = ww;
         }    
      }    
   }

   double w1   = 0.0; 
//   double w2   = 0.0; 
   double r1   = 0.0;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         int jk = j+Nt*k;
         for( i=1 ; i<Nr[jk]-1 ; ++i ){
            double ww = get_vr( theCells[jk][i].prim);
            double rr = theCells[jk][i].riph - theCells[jk][i].dr;
            //double crash_rate = ww/(r_outer-rr);
            //if( w1 < ww ) w1 = ww;
            //if( w2 < crash_rate ){ w2=crash_rate; r1 = rr; }
            if( ww > .01*max_vel && r1<rr ){ w1 = ww; r1 = rr; }
         }    
      }    
   }
   w1 = max_vel;
   struct { double value; int index; } maxw;
   maxw.value = w1;
   maxw.index = theDomain->rank;
   MPI_Allreduce( MPI_IN_PLACE , &maxw , 1 , MPI_DOUBLE_INT , MPI_MAXLOC , theDomain->theComm );
   MPI_Bcast( &r1 ,  1 , MPI_DOUBLE , maxw.index , theDomain->theComm );

   w1 = maxw.value;

   double rmin = theCells[0][0      ].riph;
   double rmax = theCells[0][Nr[0]-1].riph;
   double rhalf = ShockPos*rmax + (1.-ShockPos)*rmin;
   w1 *= ( 1. + log(r1/rhalf)/log(rmax/rhalf) );

   if( homologous_BCs ){
      w1 = 1.0;
      r1 = w1*theDomain->t;
   }
   
   if( w1 > 0.0 ){


      double w_out = w1*(rmax/r1);
      double w_in  = w1*(rmin/r1);
//This shit should not be hard-coded.  Fix this when you fix that Nickel crap. 
//      if( w_in > 1e-3 ) w_in = 1e-3;

      double rmin_new = rmin + w_in*dt;
      double rmax_new = rmax + w_out*dt;
      for( j=0 ; j<Nt ; ++j ){
         for( k=0 ; k<Np ; ++k ){

      //
      //Here is where you move the goalposts,
      //delete zones from the inside,

            int jk = j+Nt*k;
            int i=1;
            while( i<Nr[jk]-1 && theCells[jk][i].riph < rmin_new ) ++i;
            int isplit = i;
            //Add conserved stuff from inner zones into zone 0
            for( i=1 ; i<isplit ; ++i ){
               int q;
               for( q=0 ; q<NUM_Q ; ++q ){
                  theCells[jk][0].cons[q]   += theCells[jk][i].cons[q];
                  theCells[jk][0].RKcons[q] += theCells[jk][i].RKcons[q];
               }
            }
            double rm = theCells[jk][isplit-1].riph;
            double rp = theCells[jk][isplit].riph;
            double xp[3] = {rmin_new,t_jph[j]  ,p_kph[k]  };
            double xm[3] = {rm      ,t_jph[j-1],p_kph[k-1]};
            double dVl = get_dV(xp,xm);
            xp[0] = rp;
            xm[0] = rmin_new;
            double dVr = get_dV(xp,xm);
            double fracL = dVl/(dVl+dVr);
            double fracR = dVr/(dVl+dVr);
            int q;
            for( q=0 ; q<NUM_Q ; ++q ){
               theCells[jk][0].cons[q]   += fracL*theCells[jk][isplit].cons[q];
               theCells[jk][0].RKcons[q] += fracL*theCells[jk][isplit].cons[q];
            }
            if( rp-rmin_new < rmin_new-rm ){
               //if( jk==0 && theDomain->rank==0 ) printf("Killed\n");
               //if drl>drr, take the right side and add it to zone 2.
               for( q=0 ; q<NUM_Q ; ++q ){
                  theCells[jk][isplit+1].cons[q]   += fracR*theCells[jk][isplit].cons[q];
                  theCells[jk][isplit+1].RKcons[q] += fracR*theCells[jk][isplit].RKcons[q];
               }
               //In this case, you'll be deleting isplit as well.
               ++isplit;
            }else{
               //Otherwise, just account for the loss in isplit.
               for( q=0 ; q<NUM_Q ; ++q ){
                  theCells[jk][isplit].cons[q]   *= fracR;
                  theCells[jk][isplit].RKcons[q] *= fracR;
               }
            }
            int lostzones = isplit-1;
            int blocksize = Nr[jk]-isplit;
            memmove( theCells[jk]+1 , theCells[jk]+isplit , blocksize*sizeof(struct cell) );
            struct cell * c0 = theCells[jk]+0;
            struct cell * c1 = theCells[jk]+1;
            c0->riph = rmin_new;
            Nr[jk] -= lostzones;
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Nr[jk]*sizeof(struct cell) );
            c0 = theCells[jk]+0;
            c1 = theCells[jk]+1;

            xm[0] = 0.0;
            xp[0] = rmin_new;
            double rl = (2./3.)*rmin_new;
            dVl = get_dV(xp,xm);
            cons2prim( c0->cons , c0->prim , rl , .5*(xp[1]+xm[1]) , dVl );

            xm[0] = rmin_new;
            xp[0] = c1->riph;
            rp = xp[0];
            rm = xm[0];
            double rr = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
            dVr = get_dV(xp,xm);
            cons2prim( c1->cons , c1->prim , rr , .5*(xp[1]+xm[1]) , dVr );
            c1->dr = rp-rm;

      //add zones to the outside.
            theCells[jk][Nr[jk]-1].riph = rmax_new;
            double r_prev = theCells[jk][Nr[jk]-2].riph;
            double dr = rmax_new - r_prev;
            theCells[jk][Nr[jk]-1].dr = dr;
            double dtheta = t_jph[j]-t_jph[j-1];
            if( Nt==1 ) dtheta /= M_PI*theDomain->theParList.Num_R;
            int dN = (int)( (dr/rmax_new/dtheta) - 1.0 );//.5*(double)theDomain->theParList.Num_R * dr/rmax_new/log(rmax_new/rmin_new) );
            if( dN<0 ) dN = 0;
            int Nold = Nr[jk];
            Nr[jk] += dN;
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Nr[jk]*sizeof(struct cell) );
            for( i=0 ; i<dN+1 ; ++i ){
               int n = Nold+i-1;
               struct cell * c = theCells[jk]+n;
               clear_cell( c );
               double rp = r_prev + ((double)(i+1))*dr/((double)(dN+1));
               double rm = r_prev + ((double)(i  ))*dr/((double)(dN+1));
               c->riph = rp;
               c->dr   = rp-rm;
               double xp[3] = {rp,t_jph[j]  ,p_kph[k]  };
               double xm[3] = {rm,t_jph[j-1],p_kph[k-1]};
               double x[3];
               int d;
               for( d=0 ; d<3 ; ++d ) x[d] = .5*(xp[d]+xm[d]);
               initial( c->prim , x ); 
               double dV = get_dV(xp,xm);
               double rr = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
               prim2cons( c->prim , c->cons , rr , .5*(xp[1]+xm[1]) , dV );
            }
         } 
      }
   }
}

void regrid( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double jthresh = 0.1;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
   for( k=0 ; k<Np ; ++k ){
      int jk = j+Nt*k;
      for( i=2 ; i<Nr[jk]-1 ; ++i ){
         double jump = .2*log( theCells[jk][i-1].prim[RHO]/theCells[jk][i+1].prim[RHO] );
         if( fabs(jump) > jthresh ){
            struct cell * c = theCells[jk]+i;
            double rp = c->riph;
            double rm = (c-1)->riph;
            int Nsplit = (int)(fabs(jump/jthresh));
            if(Nsplit>3) Nsplit=3;
            if(Nsplit>3) printf("r=%.2e jump=%.2e split = %d (No Cause For Alarm)\n",theCells[jk][i].riph,jump,Nsplit);
            int blocksize = (Nr[jk]-1) - i;
            Nr[jk] += Nsplit;
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Nr[jk]*sizeof(struct cell) );
            memmove( theCells[jk]+i+1+Nsplit , theCells[jk]+i+1 , blocksize*sizeof(struct cell) );
            int l;
            for( l=0 ; l<Nsplit+1 ; ++l ){
               int m = l+i;
               struct cell * cnew = theCells[jk]+m;
               clear_cell( cnew );
               double rlp = rm + (rp-rm)*( (double)(l+1)/(double)(Nsplit+1) );
               double rlm = rm + (rp-rm)*( (double)(l  )/(double)(Nsplit+1) );
               double xp[3] = {rlp,t_jph[j]  ,p_kph[k]  };
               double xm[3] = {rlm,t_jph[j-1],p_kph[k-1]};
               double x[3];
               int d;
               for( d=0 ; d<3 ; ++d ) x[d] = .5*(xp[d]+xm[d]);
               initial( cnew->prim , x );
               cnew->riph = rlp;
               cnew->dr = rlp-rlm;
               double dV = get_dV(xp,xm);
               double rr = (2./3.)*(rlp*rlp*rlp-rlm*rlm*rlm)/(rlp*rlp-rlm*rlm);
               prim2cons( cnew->prim , cnew->cons , rr , x[1] , dV );
            }
            i += Nsplit;
         }
      }
   }
   }
}

void adjust_RK_cons( struct domain * theDomain , double RK ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   int i,j,k,q;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         for( i=0 ; i<Nr[j+Nt*k] ; ++i ){
            struct cell * c = &(theCells[j+Nt*k][i]);
            for( q=0 ; q<NUM_Q ; ++q ){
               c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
            }
         }
      }
   }
}

void move_cells( struct domain * theDomain , double RK , double dt){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         for( i=0 ; i<Nr[j+Nt*k] ; ++i ){
            struct cell * c = &(theCells[j+Nt*k][i]);
            c->riph = (1.-RK)*c->riph + RK*c->RKriph + c->wiph*dt;
         }
      }
   }
}

void calc_dr( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         for( i=0 ; i<Nr[j+Nt*k]-1 ; ++i ){
            double rm = 0.0; 
            if( i>0 ) rm = theCells[j+Nt*k][i-1].riph;
            double rp = theCells[j+Nt*k][i].riph;
            theCells[j+Nt*k][i].dr = rp-rm;
         }
      }
   }
}

void calc_prim( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         for( i=0 ; i<Nr[j+Nt*k] ; ++i ){
            struct cell * c = &(theCells[j+Nt*k][i]);
            double rm = 0.0; 
            if( i>0 ) rm = theCells[j+Nt*k][i-1].riph;
            double rp = c->riph;
            double r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
            double xp[3] = {rp,t_jph[j]  ,p_kph[k]  };
            double xm[3] = {rm,t_jph[j-1],p_kph[k-1]};
            double dV = get_dV( xp , xm );
            cons2prim( c->cons , c->prim , r , .5*(xp[1]+xm[1]) , dV );
         }
      }
   }
}

void plm_r( struct domain * );
void riemann_r( struct cell * , struct cell * , double , double , double );

void radial_flux( struct domain * theDomain , double dt ){
   struct cell ** theCells = theDomain->theCells;
   int Np = theDomain->Np;
   int Nt = theDomain->Nt;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int i,j,k;
   plm_r( theDomain );
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         for( i=0 ; i<Nr[j+Nt*k]-1 ; ++i ){
            struct cell * cp = theCells[j+Nt*k];
            double rp = cp[i].riph;
            double xp[3] = {rp,t_jph[j]  ,p_kph[k]  };
            double xm[3] = {rp,t_jph[j-1],p_kph[k-1]};
            double dA = get_dA(xp,xm,0); 
            riemann_r( &(cp[i]) , &(cp[i+1]) , cp[i].riph , .5*(xp[1]+xm[1]) , dA*dt );
         }
      }
   }
}

void buildfaces( struct domain * , struct face * , int * , int , int );
void plm_trans( struct domain * , struct face * , int , int );
void riemann_trans( struct face * , double , int );

void trans_flux( struct domain * theDomain , struct face * theFaces , int Nf , double dt , int dim ){
   plm_trans( theDomain , theFaces , Nf , dim );
   int f;
   for( f=0 ; f<Nf ; ++f ){
      riemann_trans( &(theFaces[f]) , dt , dim );
   }
}

int get_num_tpFaces( int , int , int );

void setup_faces( struct domain * theDomain , struct face ** theFaces , int * nn , int dim ){

   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int NN = get_num_tpFaces( Nt , Np , dim );

   buildfaces( theDomain , NULL      , nn , dim , 0 );
   int Nf = nn[NN];
   *theFaces = (struct face *) malloc( Nf*sizeof(struct face) );
   buildfaces( theDomain , *theFaces , nn , dim , 1 );

}

void source( double * , double * , double * , double * , double );

void add_source( struct domain * theDomain , double dt ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         for( i=0 ; i<Nr[j+Nt*k] ; ++i ){
            struct cell * c = &(theCells[j+Nt*k][i]);
            double rp = c->riph;
            double rm = c->riph - c->dr;
            double xp[3] = {rp,t_jph[j]  ,p_kph[k]  };
            double xm[3] = {rm,t_jph[j-1],p_kph[k-1]};
            double dV = get_dV(xp,xm);
            source( c->prim , c->cons , xp , xm , dV*dt );
         }    
      }    
   }   
}

double get_maxv( double * , double , int );

void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS , struct cell * sweep , int j , int k ){
   int Nt = theDomain->Nt;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   int jk = j + Nt*k;
   double dt = t_jph[j]-t_jph[j-1];
   *L = 0.0; 
   *S = 0.0; 
   int i;
   for( i=1 ; i<Nr[jk]-1 ; ++i ){
      struct cell * c = sweep+i;
      //double * prim = c->prim;
      //double w = c->wiph;
      double dx = c->riph*dt;///get_maxv(prim,w,1);
      double dy = c->dr;///get_maxv(prim,w,0);
      double l = dy/dx;
      double s = dx/dy;
      if( *L < l ){ *L = l; *iL = i; } 
      if( *S < s ){ *S = s; *iS = i; } 
   }
   if( Nt == 1 ){
      double rmax = sweep[Nr[jk]-1].riph;
      double rmin = sweep[0].riph;
      *S /= theDomain->theParList.Num_R*M_PI/log(rmax/rmin);
      *L *= theDomain->theParList.Num_R*M_PI/log(rmax/rmin);
   }
}

void AMRsweep( struct domain * theDomain , struct cell ** swptr , int j , int k ){
   int Nt = theDomain->Nt;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int * Nr = theDomain->Nr;

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * sweep = *swptr;
   int jk = j+Nt*k;
   double L,S;
   int iL=0;
   int iS=0;
   longandshort( theDomain , &L , &S , &iL , &iS , sweep , j , k );

   if( S>MaxShort ){

      //Possibly shift iS backwards by 1
      double drL = sweep[iS-1].dr;
      double drR = sweep[iS+1].dr;
      if( (drL < drR && iS!=1) ) --iS;

      //Remove Zone at iS+1
      sweep[iS].dr += sweep[iS+1].dr;
      sweep[iS].riph = sweep[iS+1].riph;
      sweep[iS].riph = sweep[iS+1].RKriph;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iS].cons[q]   += sweep[iS+1].cons[q];
         sweep[iS].RKcons[q] += sweep[iS+1].RKcons[q];
      }
      double rm = sweep[iS-1].riph;
      double rp = sweep[iS+1].riph;
      double r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
      double xp[3] = {rp,t_jph[j]  ,p_kph[k]  };
      double xm[3] = {rm,t_jph[j-1],p_kph[k-1]};
      double dV = get_dV( xp , xm );
      cons2prim( sweep[iS].cons , sweep[iS].prim , r , .5*(xp[1]+xm[1]) , dV );
      //Shift Memory
      int blocksize = Nr[jk]-iS-2;
      memmove( sweep+iS+1 , sweep+iS+2 , blocksize*sizeof(struct cell) );
      Nr[jk] -= 1;
      *swptr = (struct cell *) realloc( sweep , Nr[jk]*sizeof(struct cell) );
      sweep = *swptr;
      if( iS < iL ) iL -= 1;

   }

   if( L>MaxLong ){
      Nr[jk] += 1;
      *swptr = (struct cell *) realloc( sweep , Nr[jk]*sizeof(struct cell) );
      sweep = *swptr;
      int blocksize = Nr[jk]-iL-1;
      memmove( sweep+iL+1 , sweep+iL , blocksize*sizeof(struct cell) );

      double rp = sweep[iL].riph;
      double rm = sweep[iL-1].riph;
      double r0 = pow( .5*(rp*rp*rp + rm*rm*rm) , 1./3. );//sqrt( .5*(rp*rp+rm*rm) );

      sweep[iL].riph   = r0;
      sweep[iL].RKriph = r0;
      sweep[iL].dr     = r0-rm;
      sweep[iL+1].dr   = rp-r0;

      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iL].cons[q]     *= .5;
         sweep[iL].RKcons[q]   *= .5;
         sweep[iL+1].cons[q]   *= .5;
         sweep[iL+1].RKcons[q] *= .5;
      }

      rm = sweep[iL-1].riph;
      rp = sweep[iL].riph;
      double r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
      double xp[3] = {rp,t_jph[j]  ,p_kph[k]  };
      double xm[3] = {rm,t_jph[j-1],p_kph[k-1]};
      double dV = get_dV( xp , xm );
      cons2prim( sweep[iL].cons , sweep[iL].prim , r , .5*(xp[1]+xm[1]) , dV );

      rm = sweep[iL].riph;
      rp = sweep[iL+1].riph;
      r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
      xp[0] = rp;
      xm[0] = rm;
      dV = get_dV( xp , xm );
      cons2prim( sweep[iL+1].cons , sweep[iL+1].prim , r , .5*(xp[1]+xm[1]) , dV );

   }

}

void AMR( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int j,k;
   for( j=0 ; j<Nt ; ++j ){
      for( k=0 ; k<Np ; ++k ){
         AMRsweep( theDomain , theCells+j+Nt*k , j , k );
      }
   }
}

double get_entropy( double * );

void reset_entropy( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;

   if( NUM_N>0 ){
      int i,jk;
      for( jk=0 ; jk<Nt*Np ; ++jk ){
         for( i=0 ; i<Nr[jk]-1 ; ++i ){
            struct cell * c = &(theCells[jk][i]);  
            double s = get_entropy( c->prim );
            c->prim[NUM_C] = s;
            c->cons[NUM_C] = s*c->cons[DEN];
         }
      }    
   }   
}

void make_nickel( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;

   if( NUM_N>0 ){
      int i,jk;
      for( jk=0 ; jk<Nt*Np ; ++jk ){
         for( i=0 ; i<Nr[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);  
            double Pp = c->prim[PPP];
            if( c->prim[NUM_C] < Pp ){
               c->prim[NUM_C] = Pp;
               c->cons[NUM_C] = c->prim[NUM_C]*c->cons[DEN];
            }
/*
            if( Pp > 60. ){
               c->prim[NUM_C] = 1.0; 
               c->cons[NUM_C] = c->cons[DEN];
            }
*/
         }    
      }    
   }  
}

