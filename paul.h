enum{RHO,PPP,UU1,UU2};
enum{DEN,TAU,SS1,SS2};
enum{C_FIXED,C_WRIEMANN,C_WCELL};

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MOVE_CELLS C_WCELL

#define NUM_C 4
#define NUM_N 1
#define NUM_Q (NUM_C+NUM_N)
#define NUM_G 2

struct param_list{
   double t_min, t_max;
   int Num_R, Num_T, Num_P;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   double rmin, rmax;
   double thmin, thmax;
   double phimax;

   double LogZoning, MaxShort, MaxLong, Target_X, Target_W, ShockPos;
   int Absorb_BC,Move_BCs,Initial_Regrid,Initial_Cons;
   int Reset_Entropy;

   double CFL, PLM;
   double Density_Floor, Pressure_Floor;

   double Adiabatic_Index;
   int Gravity_Switch, Output_Mass;
   double PointMass;

   int Nozzle_Switch;
   double Nozzle_Power, Nozzle_Gamma, Nozzle_Eta, Nozzle_r0, Nozzle_th0, Nozzle_Time;

   double Explosion_Energy;

   int restart_flag;
};

struct domain{
   struct cell ** theCells;
   int * Nr;
   int Nt,Np,Ng;
   double * t_jph;
   double * p_kph;

   double g_point_mass;

   int rank,size;
   int dim_rank[2];
   int dim_size[2];
   int left_rank[2];
   int right_rank[2];
   MPI_Comm theComm;

   struct param_list theParList;

   double t;
   double t_init, t_fin;
   int nrpt;
   int N_rpt;
   int nsnp;
   int N_snp;
   int nchk;
   int N_chk;

   int final_step;
};

struct cell{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double grad[NUM_Q];
   double gradr[NUM_Q];
   double riph;
   double RKriph;
   double dr;
   double wiph;
};

struct face{
   struct cell * L;
   struct cell * R;
   double dxL;
   double dxR;
   double cm[3];
   double dr;
   double dA;
};

