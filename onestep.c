
#include "paul.h"

void AMR( struct domain * ); 
void move_BCs( struct domain * , double );

void clear_w( struct domain * );
void set_wcell( struct domain * );
//void set_w0( struct domain * );

void adjust_RK_cons( struct domain * , double );
void move_cells( struct domain * , double , double );
void calc_dr( struct domain * );
void calc_prim( struct domain * );

void setup_faces( struct domain * , struct face ** , int * , int );
void radial_flux( struct domain * , double dt );
void trans_flux( struct domain * , struct face * , int , double dt , int );
void add_source( struct domain * , double dt );

void boundary_r( struct domain * );
void boundary_trans( struct domain * , struct face * , int * , int );
void exchangeData( struct domain * , int );

void gravity_setup( struct domain * );
void gravity_addsrc( struct domain * , double );
void nozzle( struct domain * , double );
void add_cooling( struct domain * , double );
void make_nickel( struct domain * );

int get_num_tpFaces( int , int , int );

void onestep( struct domain * theDomain , double RK , double dt , int last_step , double global_dt ){

   int Nt = theDomain->Nt;
   int Np = theDomain->Np;

   clear_w( theDomain );
   if( MOVE_CELLS == C_WCELL ) set_wcell( theDomain );
   //set_w0( theDomain );
   int grav_switch = theDomain->theParList.Gravity_Switch;
   int add_cool    = theDomain->theParList.Add_Cooling;
   int make_nick   = theDomain->theParList.Make_Nickel;
   int nozzle_switch = theDomain->theParList.Nozzle_Switch;
   if( grav_switch ) gravity_setup( theDomain );
   adjust_RK_cons( theDomain , RK );

   struct face * theFaces_1 = NULL;
   struct face * theFaces_2 = NULL;

   int NTP1 = get_num_tpFaces( Nt , Np , 1 );
   int NTP2 = get_num_tpFaces( Nt , Np , 2 );
   int nft[ NTP1+1 ];//(Nt-1)*Np+1 ];
   int nfp[ NTP2+1 ];//Nt*(Np+1)+1 ];
   radial_flux( theDomain , dt );

   if( Nt > 1 ){
      setup_faces( theDomain , &theFaces_1 , nft , 1 );
      trans_flux( theDomain , theFaces_1 , nft[NTP1] , dt , 1 );
   }
   if( Np > 1 ){
      setup_faces( theDomain , &theFaces_2 , nfp , 2 );
      trans_flux( theDomain , theFaces_2 , nfp[NTP2] , dt , 2 );
   }
   add_source( theDomain , dt );
   if( grav_switch ) gravity_addsrc( theDomain , dt );
   if( nozzle_switch ) nozzle( theDomain , dt );
   if( add_cool ) add_cooling( theDomain , dt );
   if( make_nick ) make_nickel( theDomain );
   move_cells( theDomain , RK , dt );
   calc_dr( theDomain );
   calc_prim( theDomain );

   if( last_step ){
      AMR( theDomain );
      if( theDomain->theParList.Move_BCs ) move_BCs( theDomain , global_dt );
   }
   boundary_r( theDomain );

   if( Nt > 1 ){
      exchangeData( theDomain , 0 ); //Reverse these two when you get the chance!
      boundary_trans( theDomain , theFaces_1 , nft , 1 );
   }
   if( Np > 1 ){
      exchangeData( theDomain , 1 );
//      boundary_trans( theDomain , theFaces_2 , nfp , 2 );
   }
   if( theFaces_1 ) free( theFaces_1 );
   if( theFaces_2 ) free( theFaces_2 );
}

