
#include "../paul.h"
#include <hdf5.h>

void createFile( char * fname ){
   hid_t h5file = H5Fcreate( fname , H5F_ACC_TRUNC , H5P_DEFAULT , H5P_DEFAULT );
   H5Fclose( h5file );
}

void createGroup( char * fname , char * gname ){
   hid_t h5file = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5group = H5Gcreate1( h5file , gname , 0 );
   H5Gclose( h5group );
   H5Fclose( h5file );
}

void createDataset( char * fname , char * gname , char * dname , int dim , hsize_t * fdims , hid_t type ){
   hid_t h5file  = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5group = H5Gopen1( h5file , gname );
   hid_t fspace  = H5Screate_simple(dim,fdims,NULL);
   hid_t h5dset  = H5Dcreate1( h5group , dname , type , fspace , H5P_DEFAULT );
   H5Sclose( fspace );
   H5Dclose( h5dset );
   H5Gclose( h5group );
   H5Fclose( h5file );
}

void writeSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dwrite( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void writePatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
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

   H5Dwrite( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

int Cell2Doub( struct cell * c , double * Q , int mode ){
   if( mode==0 ) return(NUM_Q+1); else{
      int q;
      for( q=0 ; q<NUM_Q ; ++q ) Q[q] = c->prim[q];
      Q[NUM_Q] = c->riph;
      return(0);
   }
}

double get_dV( double * , double * );
void prim2cons( double * , double * , double , double , double );
void cons2prim( double * , double * , double , double , double );
void reset_entropy( struct domain * );

void output( struct domain * theDomain , char * filestart ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   int Nt_Tot = theDomain->theParList.Num_T;
   int Np_Tot = theDomain->theParList.Num_P;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;


   char filename[256];
   sprintf(filename,"%s.h5",filestart);
//   char filename_1d[256];
//   sprintf(filename_1d,"%s_1d.dat",filestart);
   int jmin = 0;
   if( dim_rank[0] != 0 ) jmin = Ng;
   int jmax = Nt;
   if( dim_rank[0] != dim_size[0]-1 ) jmax = Nt-Ng;
/*
   int kmin = 0;
   int kmax = Np;
   if( dim_rank[1] != 0 ) kmin = Ng;
   if( dim_rank[1] != dim_size[1]-1 ) kmax = Np-Ng;
*/
   int kmin = Ng;
   int kmax = Np-Ng;
   if( Np == 1 ){ kmin = 0 ; kmax = 1; }

   int Ntot = 0;
   int j,k;
   for( j=jmin ; j<jmax ; ++j ){
      for( k=kmin ; k<kmax ; ++k ){
         int jk = j+Nt*k;
         Ntot += Nr[jk];
      }
   }
   int myNtot = Ntot;
   MPI_Allreduce( MPI_IN_PLACE , &Ntot , 1 , MPI_INT , MPI_SUM , theDomain->theComm );

   int Ndoub = Cell2Doub(NULL,NULL,0);

   hsize_t fdims1[1];
   hsize_t fdims2[2];
   if( rank==0 ){
      printf("Writing Checkpoint...\n");
      
      createFile(filename);
      createGroup(filename,"Grid");

      fdims1[0] = 1;
      createDataset(filename,"Grid","T",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims1[0] = Nt_Tot+1;
      createDataset(filename,"Grid","t_jph",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims1[0] = Np_Tot+1;
      createDataset(filename,"Grid","p_kph",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims2[0] = Nt_Tot;
      fdims2[1] = Np_Tot;
      createDataset(filename,"Grid","Index",2,fdims2,H5T_NATIVE_INT);
      createDataset(filename,"Grid","Nr",2,fdims2,H5T_NATIVE_INT);

      createGroup(filename,"Data");

      fdims2[0] = Ntot;
      fdims2[1] = Ndoub;
      createDataset(filename,"Data","Cells",2,fdims2,H5T_NATIVE_DOUBLE);
   }
   MPI_Barrier( theDomain->theComm );
   if( rank==0 ){
      writeSimple(filename,"Grid","T",&(theDomain->t),H5T_NATIVE_DOUBLE);
   }

   int nrk;
   int j0 = 0;
   int jSum = 0;
   int k0 = 0;
   int kSum = 0;
/*
   for( nrk=0 ; nrk < size ; ++nrk ){
      if( nrk == rank ){
         j0 = jSum;
         k0 = kSum;
         jSum += jmax-jmin;
         kSum += kmax-kmin;
      }else{ jSum=0 ; kSum = 0; }
      MPI_Allreduce( MPI_IN_PLACE , &jSum , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
      MPI_Allreduce( MPI_IN_PLACE , &kSum , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   }
*/
   for( nrk=0 ; nrk < dim_size[0] ; ++nrk ){
      if( nrk == dim_rank[0] ){
         j0 = jSum;
         if( dim_rank[1] == 0 ) jSum += jmax-jmin; else jSum = 0;
      }else{jSum=0;}
      MPI_Allreduce( MPI_IN_PLACE , &jSum , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   }
   for( nrk=0 ; nrk < dim_size[1] ; ++nrk ){
      if( nrk == dim_rank[1] ){
         k0 = kSum;
         if( dim_rank[0] == 0 ) kSum += kmax-kmin; else kSum = 0; 
      }else{kSum=0;}
      MPI_Allreduce( MPI_IN_PLACE , &kSum , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   }
//printf("rank=%d,j0=%d,jSum=%d, k0=%d,kSum=%d\n",rank,j0,jSum,k0,kSum);
   if( Nt_Tot == 1 ){ j0 = 0; jSum = 1; }
   if( Np_Tot == 1 ){ k0 = 0; kSum = 1; }
   int myIndex=0;
   for( nrk=0 ; nrk < size ; ++nrk ){
   if( rank==nrk ){
      int jSize = jmax-jmin;
      int kSize = kmax-kmin;
      int Index[jSize*kSize];
      int Size[jSize*kSize];
      double * Qwrite = (double *) malloc( myNtot*Ndoub*sizeof(double) );
      int index = 0;
      for( k=kmin ; k<kmax ; ++k ){
         for( j=jmin ; j<jmax ; ++j ){
            //int jk = (j-jmin) + jSize*(k-kmin);
            int jk = (j-jmin)*kSize + (k-kmin);
            Index[jk] = index + myIndex;
            Size[jk]  = Nr[j+Nt*k];
            int i;
            for( i=0 ; i<Nr[j+Nt*k] ; ++i ){
               struct cell * c = &(theCells[j+Nt*k][i]);
               Cell2Doub( c , Qwrite+index*Ndoub , 1 );
               index++;
            }
         }
      }
      //Write Cell Data
      int start2[2]    = {myIndex,0};
      int loc_size2[2] = {myNtot,Ndoub};
      int glo_size2[2] = {Ntot,Ndoub};
      writePatch( filename , "Data" , "Cells" , Qwrite , H5T_NATIVE_DOUBLE , 2 , start2 , loc_size2 , glo_size2 );
      free(Qwrite);
      //Write Indices and Sizes for each radial track
      start2[0] = j0;
      start2[1] = k0;
      loc_size2[0] = jSize;
      loc_size2[1] = kSize;
      glo_size2[0] = Nt_Tot;
      glo_size2[1] = Np_Tot;
      writePatch( filename , "Grid" , "Index" , Index , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
      writePatch( filename , "Grid" , "Nr"    , Size  , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
      //Write 1D Theta Data
      if( dim_rank[1] == 0 ){
         int offset = Ng;
         if( dim_rank[0] == 0 ) offset = 0;
         int start1[1]    = {j0};
         int loc_size1[1] = {jSize};
         if( dim_rank[0] == dim_size[0]-1 ) loc_size1[0]++;
         int glo_size1[1] = {Nt_Tot+1};
         writePatch( filename , "Grid" , "t_jph" , t_jph-1+offset , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
      }
      //Write 1D Phi Data
      if( dim_rank[0] == 0 ){
         int offset = Ng;
         if( dim_rank[1] == 0 ) offset = 0;
         int start1[1]    = {k0};
         int loc_size1[1] = {kSize};
         if( dim_rank[1] == dim_size[1]-1 ) loc_size1[0]++;
         int glo_size1[1] = {Np_Tot+1};
         writePatch( filename , "Grid" , "p_kph" , p_kph-1+offset , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
      }
   }
      myIndex += myNtot;
      MPI_Bcast( &myIndex , 1 , MPI_INT , nrk , theDomain->theComm );
   }

   if( theDomain->theParList.Reset_Entropy ) reset_entropy( theDomain );
/*
   int j_min = 0;
   int j_max = Nt;
   int k_min = 0;
   int k_max = Np;

   if( dim_rank[0] != 0 ) j_min = NUM_G;
   if( dim_rank[0] != dim_size[0]-1 ) j_max = Nt-NUM_G;
   if( dim_rank[1] != 0 ) k_min = NUM_G;
   if( dim_rank[1] != dim_size[1]-1 ) k_max = Np-NUM_G;

   int rk;
   for( rk=0 ; rk<size ; ++rk ){
      if( rank==rk ){
         FILE * pFile = fopen( filename , "a" );
         int i,j,k,q;
         for( j=j_min ; j<j_max ; ++j ){
            double th  = .5*(t_jph[j]+t_jph[j-1]);
            double dth = t_jph[j]-t_jph[j-1];
            for( k=k_min ; k<k_max ; ++k ){
               double phi = .5*(p_kph[k]+p_kph[k-1]);
               double dph = p_kph[k]-p_kph[k-1];
               int jk = j+Nt*k;
               for( i=0 ; i<Nr[jk] ; ++i ){
                  struct cell * c = &(theCells[jk][i]);
                  double r  = c->riph - .5*c->dr;
                  double dr = c->dr;
if( rank==0 && i==0 ) fprintf(pFile,"# ");
                  fprintf(pFile,"%e %e %e ",r,th,phi);
                  fprintf(pFile,"%e %e %e ",dr,dth,dph);
                  for( q=0 ; q<NUM_Q ; ++q ){
                     fprintf(pFile,"%e ",c->prim[q]);
                  }
                  fprintf(pFile,"\n");
               }
            }
         }
         fclose( pFile );
      }
      MPI_Barrier( theDomain->theComm );
   }

   double Rmin = theCells[0][0].riph;
   double Rmax = theCells[0][Nr[0]-1].riph;
 
   int NUM_R = theDomain->theParList.Num_R;
   int i,j,k,q;
   double cons_1d_avg[NUM_R*NUM_Q];
   double prim_1d_avg[NUM_R*NUM_Q];
   for( i=0 ; i<NUM_R*NUM_Q ; ++i ){
      cons_1d_avg[i] = 0.0;
      prim_1d_avg[i] = 0.0;
   }
   
   for( j=j_min ; j<j_max ; ++j ){
      double thp = t_jph[j];
      double thm = t_jph[j-1];
      for( k=k_min ; k<k_max ; ++k ){
         double php = p_kph[k];
         double phm = p_kph[k-1];
         int jk = j+Nt*k;
         for( i=1 ; i<Nr[jk] ; ++i ){
            double rp = theCells[jk][i].riph;
            double rm = rp-theCells[jk][i].dr;
            double ip = (NUM_R)*(rp-Rmin)/(Rmax-Rmin);//(NUM_R)*log(rp/Rmin)/log(Rmax/Rmin);
            double im = (NUM_R)*(rm-Rmin)/(Rmax-Rmin);//(NUM_R)*log(rm/Rmin)/log(Rmax/Rmin);
            double xp[3] = {rp,thp,php};
            double xm[3] = {rm,thm,phm};
            double dV = get_dV( xp , xm );
          
            double delta_i = ip-im; 
            int iip = (int) ip;
            int iim = ((int) im);
            double dip = ip-(double)iip;
            double dim = (double)(iim+1)-im;

            if( iip == iim ) dip = delta_i;
        
            int iCons;
            for( iCons = iim ; iCons <= iip ; ++iCons ){
               double frac = 1./delta_i;
               if( iCons == iim ) frac = dim/delta_i;
               if( iCons == iip ) frac = dip/delta_i;

               int i2 = iCons;
               if( i2 < 0 ) i2 = 0;
               if( i2 > NUM_R-1 ) i2 = NUM_R-1;

               for( q=0 ; q<NUM_Q ; ++q ){
                  cons_1d_avg[i2*NUM_Q + q] += frac*theCells[jk][i].cons[q];
                  prim_1d_avg[i2*NUM_Q + q] += frac*theCells[jk][i].prim[q]*dV;
               }
            }
         }
      }
   }

   MPI_Allreduce( MPI_IN_PLACE , cons_1d_avg , NUM_R*NUM_Q , MPI_DOUBLE , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , prim_1d_avg , NUM_R*NUM_Q , MPI_DOUBLE , MPI_SUM , theDomain->theComm );

   double THETA_MIN = theDomain->theParList.thmin;
   double THETA_MAX = theDomain->theParList.thmax;
   double PHI_MAX   = theDomain->theParList.phimax;

   if( rank==0 ){
      FILE * pFile_1d = fopen(filename_1d,"w");
      for( i=0 ; i<NUM_R ; ++i ){
         double P_out[NUM_Q];
         double xxm = (double)i/(double)NUM_R;
         double xxp = (double)(i+1)/(double)NUM_R;
         double rm = xxm*(Rmax-Rmin)+Rmin;//pow(Rmax/Rmin,xxm)*Rmin;
         double rp = xxp*(Rmax-Rmin)+Rmin;//pow(Rmax/Rmin,xxp)*Rmin;
         double xp[3] = {rp,THETA_MAX,PHI_MAX};
         double xm[3] = {rm,THETA_MIN,0.0};
         double r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
         double dV = get_dV( xp , xm );
         cons2prim( &(cons_1d_avg[i*NUM_Q]) , P_out , r , dV );
         fprintf(pFile_1d,"%e ",r);
         for( q=0 ; q<NUM_Q ; ++q ){
            fprintf(pFile_1d,"%e ",P_out[q]);
            fprintf(pFile_1d,"%e ",prim_1d_avg[i*NUM_Q+q]/dV);
         }
         fprintf(pFile_1d,"\n");
      }
      fclose( pFile_1d );
   }
*/
}
