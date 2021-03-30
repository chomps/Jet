
#include "../paul.h"
#include <hdf5.h>
#include <string.h>

void getH5dims( char * file , char * group , char * dset , hsize_t * dims ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );
   hid_t h5spc = H5Dget_space( h5dst );

   H5Sget_simple_extent_dims( h5spc , dims , NULL);

   H5Sclose( h5spc );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readPatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];

   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = loc_size[d];
      fdims[d] = glo_size[d];

      fstart[d]  = start[d];
      fstride[d] = 1;
      fcount[d]  = loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}


void Doub2Cell( double * Q , struct cell * c ){
   int q;
   for( q=0 ; q<NUM_Q ; ++q ) c->prim[q] = Q[q];
   c->riph = Q[NUM_Q];
}

int getN0( int , int , int );
void freeDomain( struct domain * );

void restart( struct domain * theDomain ){

   //This code has not been bug-tested in 3D.
   //The ordering of indices in particular hasn't been checked
   //i.e. should I use jk = j + Nt*k or jk = j*Np + k?

   freeDomain( theDomain );

   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   int i,j,k;

   char filename[256];
   strcpy(filename,"input.h5");
   char group1[256];
   strcpy(group1,"Grid");
   char group2[256];
   strcpy(group2,"Data");

   hsize_t dims[3];

   if( rank==0 ) printf("Restarting from file...\n");

   int NUM_T,NUM_P;
   double tstart;
   //Read the time from "T" and get Nt and Np from 
   //the dimensions of "Index".  Broadcast to all ranks.
   if(rank==0){
      readSimple( filename , group1 ,"T", &tstart , H5T_NATIVE_DOUBLE );
      getH5dims( filename , group1 ,"Index", dims );
      NUM_T = dims[0];
      NUM_P = dims[1];
   }
   MPI_Bcast( &NUM_T  , 1 , MPI_INT    , 0 , theDomain->theComm );
   MPI_Bcast( &NUM_P  , 1 , MPI_INT    , 0 , theDomain->theComm );
   MPI_Bcast( &tstart , 1 , MPI_DOUBLE , 0 , theDomain->theComm );
   theDomain->theParList.Num_T = NUM_T;
   theDomain->theParList.Num_P = NUM_P;
   theDomain->t = tstart;
 
 
   //The following is very similar to equivalent code in gridsetup.c
   //Now you're just doing the process over because you're restarting
   //from file. 
   int N0t = getN0( dim_rank[0] , dim_size[0] , NUM_T );
   int N1t = getN0( dim_rank[0]+1 , dim_size[0] , NUM_T );
   int Nt = N1t-N0t;

   int N0p = getN0( dim_rank[1] , dim_size[1] , NUM_P );
   int N1p = getN0( dim_rank[1]+1 , dim_size[1] , NUM_P );
   int Np = N1p-N0p;

   if( dim_rank[0] != 0 ){ Nt += Ng; N0t -= Ng;}
   if( dim_rank[0] != dim_size[0]-1 ) Nt += Ng;
   if( NUM_P != 1 ){ Np += 2*Ng; }//N0p -= Ng; }

   theDomain->Nt = Nt;
   theDomain->Np = Np;

   theDomain->Nr    = (int *)    malloc( Np*Nt*sizeof(int) );
   theDomain->t_jph = (double *) malloc( (Nt+1)*sizeof(double) );
   theDomain->p_kph = (double *) malloc( (Np+1)*sizeof(double) );

   theDomain->theCells = (struct cell **) malloc( Nt*Np*sizeof( struct cell * ) );

   struct cell ** theCells = theDomain->theCells;


   //The following must happen in serial because different processors
   //will try to read from the same file.  In principle we can write
   //this using parallel HDF5, but I'd rather cross that bridge 
   //when I come to it.
   int nrk;
   int Nq=0;
   for( nrk=0 ; nrk<size ; ++nrk ){
   if( rank==nrk ){

      //Read the theta values of the grid...
      int start1[1] = {N0t};
      int loc_size1[1] = {Nt+1};
      int glo_size1[1] = {NUM_T+1};
      double t_jph[Nt+1];
      readPatch( filename , group1 ,"t_jph", t_jph , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 ); 
      for( j=0 ; j<Nt+1 ; ++j ) theDomain->t_jph[j] = t_jph[j];

int kmin = 0;
int kmax = Np;
if( Np>1 ){ kmin = Ng ; kmax = Np-Ng; }
int Np_read = kmax-kmin;
      //Read the phi values of the grid...
      start1[0]    = N0p;
      loc_size1[0] = Np_read+1;
      glo_size1[0] = NUM_P+1;
      double p_kph[Np_read+1];
      readPatch( filename , group1 ,"p_kph", p_kph , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
      for( k=kmin ; k<kmax+1 ; ++k ) theDomain->p_kph[k] = p_kph[k-kmin];
if( Np>1 ){
   theDomain->p_kph[0 ] = 2.*theDomain->p_kph[   1]-theDomain->p_kph[   2];
   theDomain->p_kph[Np] = 2.*theDomain->p_kph[Np-1]-theDomain->p_kph[Np-2];
}

      ++(theDomain->t_jph);
      ++(theDomain->p_kph);

      //Read the indexing information so you know how to read in
      //The radial tracks of data which are coming up...
      int start2[2]   = {N0t,N0p};
      int loc_size2[2] = {Nt,Np_read};
      int glo_size2[2] = {NUM_T,NUM_P};
      int Nr[Nt*Np];
      int Index[Nt*Np];
      readPatch( filename , group1 ,"Nr"   , Nr    , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
      readPatch( filename , group1 ,"Index", Index , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );

      getH5dims( filename , group2 ,"Cells", dims );
      int Nc = dims[0];
      Nq = dims[1];

      //Read in each radial track one at a time, because
      //you don't really know where the different radial
      //tracks belong in memory, they might not be 
      //contiguous, and they definitely are not in a rectangular block.
      for( j=0 ; j<Nt ; ++j ){
         for( k=kmin ; k<kmax ; ++k ){
            int jk = j+Nt*k;
            int jk_read = j+Nt*(k-kmin);
            theDomain->Nr[jk] = Nr[jk_read];
            start2[0] = Index[jk_read];
            start2[1] = 0;
            loc_size2[0] = Nr[jk_read];
            loc_size2[1] = Nq;
            glo_size2[0] = Nc;
            glo_size2[1] = Nq;
            double TrackData[Nr[jk_read]*Nq];
            readPatch( filename , group2 ,"Cells", TrackData , H5T_NATIVE_DOUBLE , 2 , start2 , loc_size2 , glo_size2 );
            theDomain->theCells[jk] = (struct cell *) malloc( Nr[jk_read]*sizeof(struct cell) );
            for( i=0 ; i<Nr[jk_read] ; ++i ){
               struct cell * c = theCells[jk]+i;
               Doub2Cell( TrackData + i*Nq , c );
            }
         }
         for( k=0 ; k<kmin ; ++k ){
            int jk = j+Nt*k;
            theDomain->Nr[jk] = 3;
            theDomain->theCells[jk] = (struct cell *) malloc( 3*sizeof(struct cell) );
         }
         for( k=kmax ; k<Np ; ++k ){
            int jk = j+Nt*k;
            theDomain->Nr[jk] = 3;
            theDomain->theCells[jk] = (struct cell *) malloc( 3*sizeof(struct cell) );
         }
      }

   }
   MPI_Barrier(theDomain->theComm);
   }
   if( Nq != NUM_Q+1 ){ if(rank==0)printf("Ummm, I got an hdf5 read error. Check NUM_Q.\n"); exit(1); }

}

