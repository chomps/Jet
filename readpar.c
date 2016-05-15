enum{VAR_INT,VAR_DOUB,VAR_STR};

#include "paul.h"
#include <string.h>

int readvar( char * filename , char * varname , int vartype , void * ptr ){

   FILE * inFile = fopen( filename , "r" );
   char s[512];
   char nm[512];
   char s1[512];
   int found = 0;
   
   while( (fgets(s,512,inFile) != NULL) && found==0 ){
      sscanf(s,"%s ",nm);
      if( strcmp(nm,varname)==0 ){
         strcpy(s1,s);
         found=1;
      }
   }
   
   fclose( inFile );
   if( found==0 ){ printf("%s not found\n",varname); return(1); }

   char * s2 = s1+strlen(nm)+strspn(s1+strlen(nm),"\t :=>_");

   double temp;
   char stringval[256];

   sscanf(s2,"%lf",&temp);
   sscanf(s2,"%256s",stringval);

   if( vartype == VAR_INT ){
      *((int *)   ptr) = (int)temp;
   }else if( vartype == VAR_DOUB ){
      *((double *)ptr) = (double)temp;
   }else{
      strcpy( ptr , stringval );
   }

   return(0);
}

int read_par_file( struct domain * theDomain ){

   MPI_Comm_rank(MPI_COMM_WORLD,&(theDomain->rank));
   MPI_Comm_size(MPI_COMM_WORLD,&(theDomain->size));

   int rank = theDomain->rank;
   int size = theDomain->size;

   struct param_list * theList = &( theDomain->theParList );

   char pfile[] = "in.par";

   int err=0;  

   int nrank;
   for( nrank=0 ; nrank<size ; ++nrank ){
      if( rank==nrank ){
         err += readvar( pfile , "Num_R"           , VAR_INT  , &(theList->Num_R)           );
         err += readvar( pfile , "Num_Theta"       , VAR_INT  , &(theList->Num_T)           );
         err += readvar( pfile , "Num_Phi"         , VAR_INT  , &(theList->Num_P)           );
         err += readvar( pfile , "Num_Reports"     , VAR_INT  , &(theList->NumRepts)        );
         err += readvar( pfile , "Num_Snapshots"   , VAR_INT  , &(theList->NumSnaps)        );
         err += readvar( pfile , "Num_Checkpoints" , VAR_INT  , &(theList->NumChecks)       );
         err += readvar( pfile , "T_Start"         , VAR_DOUB , &(theList->t_min)           );
         err += readvar( pfile , "T_End"           , VAR_DOUB , &(theList->t_max)           );
         err += readvar( pfile , "R_Min"           , VAR_DOUB , &(theList->rmin)            );
         err += readvar( pfile , "R_Max"           , VAR_DOUB , &(theList->rmax)            );
         err += readvar( pfile , "Theta_Min"       , VAR_DOUB , &(theList->thmin)           );
         err += readvar( pfile , "Theta_Max"       , VAR_DOUB , &(theList->thmax)           );
         err += readvar( pfile , "Phi_Max"         , VAR_DOUB , &(theList->phimax)          );
         err += readvar( pfile , "Log_Zoning"      , VAR_DOUB , &(theList->LogZoning)       );
         err += readvar( pfile , "Max_Aspect_Short", VAR_DOUB , &(theList->MaxShort)        );
         err += readvar( pfile , "Max_Aspect_Long" , VAR_DOUB , &(theList->MaxLong)         );
         err += readvar( pfile , "CFL"             , VAR_DOUB , &(theList->CFL)             );
         err += readvar( pfile , "PLM"             , VAR_DOUB , &(theList->PLM)             );
         err += readvar( pfile , "Adiabatic_Index" , VAR_DOUB , &(theList->Adiabatic_Index) );
         err += readvar( pfile , "Density_Floor"   , VAR_DOUB , &(theList->Density_Floor)   );
         err += readvar( pfile , "Pressure_Floor"  , VAR_DOUB , &(theList->Pressure_Floor)  );
         err += readvar( pfile , "Gravity_Switch"  , VAR_INT  , &(theList->Gravity_Switch)  );
         err += readvar( pfile , "Output_Mass"     , VAR_INT  , &(theList->Output_Mass)     );
         err += readvar( pfile , "PointMass"       , VAR_DOUB , &(theList->PointMass)       );
         err += readvar( pfile , "Explosion_Energy", VAR_DOUB , &(theList->Explosion_Energy));
         err += readvar( pfile , "Gam_0"           , VAR_DOUB , &(theList->Gam_0)           );
         err += readvar( pfile , "Gam_Boost"       , VAR_DOUB , &(theList->Gam_Boost)       );
         err += readvar( pfile , "Absorbing_BC"    , VAR_INT  , &(theList->Absorb_BC)       );
         err += readvar( pfile , "Move_Boundaries" , VAR_INT  , &(theList->Move_BCs)        );
         err += readvar( pfile , "Shock_Position"  , VAR_DOUB , &(theList->ShockPos)        );
         err += readvar( pfile , "Initial_Regrid"  , VAR_INT  , &(theList->Initial_Regrid)  );
         err += readvar( pfile , "Restart"         , VAR_INT  , &(theList->restart_flag)    );
         err += readvar( pfile , "Target_X0"       , VAR_DOUB , &(theList->Target_X)        );
         err += readvar( pfile , "Target_Weight"   , VAR_DOUB , &(theList->Target_W)        );
         err += readvar( pfile , "Nozzle_Switch"   , VAR_INT  , &(theList->Nozzle_Switch)   );
         err += readvar( pfile , "Nozzle_Power"    , VAR_DOUB , &(theList->Nozzle_Power)    );
         err += readvar( pfile , "Nozzle_Gamma"    , VAR_DOUB , &(theList->Nozzle_Gamma)    );
         err += readvar( pfile , "Nozzle_Eta"      , VAR_DOUB , &(theList->Nozzle_Eta)      );
         err += readvar( pfile , "Nozzle_r0"       , VAR_DOUB , &(theList->Nozzle_r0)       );
         err += readvar( pfile , "Nozzle_th0"      , VAR_DOUB , &(theList->Nozzle_th0)      );
         err += readvar( pfile , "Nozzle_Time"     , VAR_DOUB , &(theList->Nozzle_Time)     );
         err += readvar( pfile , "Use_Logtime"     , VAR_INT  , &(theList->Out_LogTime)     );
         err += readvar( pfile , "Initial_Cons"    , VAR_INT  , &(theList->Initial_Cons)    );
         err += readvar( pfile , "Reset_Entropy"   , VAR_INT  , &(theList->Reset_Entropy)   );
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   int errtot;
   MPI_Allreduce( &err , &errtot , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );

   if( errtot > 0 ){
      printf("Read Failed, err = %d\n",err);
      return(1);
   }

   return(0);

}


