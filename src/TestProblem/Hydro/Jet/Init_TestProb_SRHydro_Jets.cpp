#include <random>
#include <limits.h>
#include <math.h>
#include "GAMER.h"
#include "TestProb.h"

void ***calloc_3d_array (size_t nt, size_t nr, size_t nc, size_t size);
void free_3d_array(void ***array);
void Mis_Cartesian2Spherical( const double Cartesian[], double Spherical[] );
void CartesianRotate( double x[], double theta, double phi, bool inverse );
void Interpolation_UM_IC( real x, real y, real z, real *Pri );
real TrilinearInterpolation(real *FieldAtVertices, real *xyz000, real *dxyz, real *xyz);
void SetArray();

#if ( MODEL == HYDRO )

// problem-specific global variables
// =======================================================================================

// options
static int      Jet_Ambient;             // [0/1/9]: uniform/Milky-Way/load-from-file
static bool     Jet_Precession;          // flag: precessing jet source
static bool     Jet_TimeDependentSrc;    // flag: time-dependent fluid variables in source
static int      Jet_Fire;                // [0/1/2/3]: no jet/jet1/jet2/bipolar jet
static double   Jet_Duration;            // a duration of jet injection from the start of simulation

// general parameters
static double   ParticleMass;            // atomic mass unit in jet source

// uniform background parameters
static double   Amb_UniformDens;         // uniform ambient density
static double   Amb_UniformVel[3];       // uniform ambient 4-velocity
static double   Amb_UniformTemp;         // uniform ambient temperature

// Milky Way parameters
       double   IsothermalSlab_Center[3];

// jet fluid parameters
static double   Jet_SrcVel;              // jet 4-velocity
static double   Jet_SrcDens;             // jet density
static double   Jet_SrcTemp;             // jet temperature
static bool     Jet_SmoothVel;           // smooth radial component of 4-velocity on cross section


// sound speed
static double   CharacteristicSpeed;     // the characteristic speed of the simulation problem
                                         // the default end-time (END_T) will be estimated from
                                         // `CharacteristicSpeed` and `BOX_SIZE`
static real *buffer;
static real *Header;
static real ***Rhoo;
static real ***VelX;
static real ***VelY;
static real ***VelZ;
static real ***Pres;
static real *X;
static real *Y;
static real *Z;

// fermi bubbles
       real IsothermalSlab_VelocityDispersion;
       real IsothermalSlab_PeakDens;
       real IsothermalSlab_Truncation;
static real ambientTemperature;
static real gasDiskTemperature;
static real gasDiskPeakDens;
static real interfaceHeight;

// Dark logarithmic halo potential
       real  v_halo;
       real  distance_h;

void Init_ExtPot_IsothermalSlab(); 

// =======================================================================================
/*        G       A       C              */ 
/*          ____________                 */
/*          \     |     /                */
/*           \   E|    /        z        */
/*            \  /|\  /         ^        */
/*             \/_|_\/          |        */
/*             /\O| /\B         |        */
/*            /  \|/  \                  */
/*           /   D|    \                 */
/*          /_____|_____\                */
/*                F                      */
// =======================================================================================
//
// jet geometry parameters
static double   Jet_Radius;              // length of OB
static double   Jet_HalfHeight;          // length of OA
static double   Jet_HalfOpeningAngle;    // half-opening angle (i.e. ∠ADC)
static double   Jet_CenOffset[3];        // jet central coordinates offset
static double   Jet_Center[3];           // jet central coordinates
static double   Jet_MaxDis;              // maximum distance between the jet source and the jet center


// precession parameters
static double   Jet_AngularVelocity;     // precession angular velocity (degree per code_time)
static double   Jet_PrecessionAngle;     // precession angle in degree
static double   Jet_PrecessionAxis[3];   // cone orientation vector (x,y,z). i.e. vector OA
                                         // --> NOT necessary to be a unit vector

// time-depent source
static double   Jet_BurstStartTime;      // start burst time in jet source
static double   Jet_BurstEndTime;        // end burst time in jet source
static double   Jet_Burst4VelRatio;      // increase 4-velocity     by a factor of `Jet_Burst4VelRatio` during `Jet_BurstStartTime` and `Jet_BurstEndTime`
static double   Jet_BurstDensRatio;      // increase proper density by a factor of `Jet_BurstDensRatio` during `Jet_BurstStartTime` and `Jet_BurstEndTime` 
static double   Jet_BurstTempRatio;      // increase temperature    by a factor of `Jet_BurstTempRatio` during `Jet_BurstStartTime` and `Jet_BurstEndTime` 
static bool     Flag_Burst4Vel;          // flag: burst 4-velocity
static bool     Flag_BurstDens;          // flag: burst proper density
static bool     Flag_BurstTemp;          // flag: burst temperature

static double   Amb_FluSphereRadius;     //

// =======================================================================================



//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifndef GRAVITY
   if ( Jet_Ambient == 1 )
        Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

// warnings
   if ( MPI_Rank == 0 )
   {
      for (int s=0; s<6; s++)
         if ( OPT__BC_FLU[s] != BC_FLU_OUTFLOW )
            Aux_Message( stderr, "WARNING : it's recommended to use the outflow BC (currently OPT__BC_FLU[%d] = %d != 2)\n",
                         s, OPT__BC_FLU[s] );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",         &VARIABLE,             DEFAULT,                   MIN,            MAX    );
// ************************************************************************************************************************

// load options
   ReadPara->Add( "Jet_Ambient",             &Jet_Ambient,              1,                       0,              9    );
   ReadPara->Add( "Jet_Fire",                &Jet_Fire,                 3,                       0,              3    );
   ReadPara->Add( "Jet_Precession",          &Jet_Precession,           false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_TimeDependentSrc",    &Jet_TimeDependentSrc,     false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_Duration",            &Jet_Duration        ,     NoMax_double,          0.0,   NoMax_double    );

// load jet fluid parameters
   ReadPara->Add( "Jet_SrcVel",              &Jet_SrcVel    ,          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_SmoothVel",           &Jet_SmoothVel ,           false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_SrcDens",             &Jet_SrcDens   ,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_SrcTemp",             &Jet_SrcTemp   ,          -1.0,          Eps_double,     NoMax_double    );

// load source geometry parameters
   ReadPara->Add( "Jet_Radius",              &Jet_Radius,              -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HalfHeight",          &Jet_HalfHeight,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HalfOpeningAngle",    &Jet_HalfOpeningAngle,    -1.0,                   0.0,           90.0    );
   ReadPara->Add( "Jet_CenOffset_x",         &Jet_CenOffset [0],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_y",         &Jet_CenOffset [1],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_z",         &Jet_CenOffset [2],        NoDef_double, NoMin_double,   NoMax_double    );

// load precission parameters
   ReadPara->Add( "Jet_AngularVelocity",     &Jet_AngularVelocity,      NoDef_double,          0.0,   NoMax_double    );
   ReadPara->Add( "Jet_PrecessionAngle",     &Jet_PrecessionAngle,      NoDef_double, NoMin_double,           90.0    );
   ReadPara->Add( "Jet_PrecessionAxis_x",    &Jet_PrecessionAxis[0],    NoDef_double, NoMin_double,   NoMax_double    );               
   ReadPara->Add( "Jet_PrecessionAxis_y",    &Jet_PrecessionAxis[1],    NoDef_double, NoMin_double,   NoMax_double    );  
   ReadPara->Add( "Jet_PrecessionAxis_z",    &Jet_PrecessionAxis[2],    NoDef_double, NoMin_double,   NoMax_double    );  

// load uniform background parameters
   ReadPara->Add( "Amb_UniformDens",         &Amb_UniformDens,         -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Amb_UniformVel_x",        &Amb_UniformVel[0],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Amb_UniformVel_y",        &Amb_UniformVel[1],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Amb_UniformVel_z",        &Amb_UniformVel[2],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Amb_UniformTemp",         &Amb_UniformTemp,         -1.0,          Eps_double,     NoMax_double    );


   ReadPara->Add( "Amb_FluSphereRadius",     &Amb_FluSphereRadius,     -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "CharacteristicSpeed",     &CharacteristicSpeed,     -1.0,          NoMin_double,   NoMax_double    );

// load Milky Way parameters
   ReadPara->Add( "IsothermalSlab_Center_x",          &IsothermalSlab_Center[0],          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "IsothermalSlab_Center_y",          &IsothermalSlab_Center[1],          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "IsothermalSlab_Center_z",          &IsothermalSlab_Center[2],          -1.0,          NoMin_double,   NoMax_double    );

// load time-dependent source varibles
   ReadPara->Add( "Jet_BurstStartTime",      &Jet_BurstStartTime,      -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstEndTime",        &Jet_BurstEndTime,        -1.0,          NoMin_double,   NoMax_double    );

   ReadPara->Read( FileName );

   delete ReadPara;

// Read header for the fermi bubbles
   SetArray();

// replace useless parameters with NaN
   if ( Jet_Ambient != 0 )
   {
     Amb_UniformDens      = NAN;
     Amb_UniformVel[0]    = NAN;
     Amb_UniformVel[1]    = NAN;
     Amb_UniformVel[2]    = NAN;
   }

   if ( Amb_FluSphereRadius < 0.0 )
   {
      Amb_FluSphereRadius = NAN;
   }

   
   if ( !Jet_TimeDependentSrc )
   {
     Jet_BurstDensRatio  = NAN;
     Jet_Burst4VelRatio  = NAN;
     Jet_BurstTempRatio  = NAN; 
     Jet_BurstStartTime  = NAN;
     Jet_BurstEndTime    = NAN;
   }

// (1-2) check runtime parameters

// check time-dependent source
   if ( Jet_TimeDependentSrc )
   {
     if ( !Flag_Burst4Vel && !Flag_BurstDens && !Flag_BurstTemp )
       Aux_Error( ERROR_INFO, "One of Flag_Burst4Vel, Flag_BurstDens or Flag_BurstTemp must be enabled !!\n" );

     if ( Jet_BurstEndTime <= Jet_BurstStartTime)
       Aux_Error( ERROR_INFO, "Jet_BurstEndTime <= Jet_BurstStartTime !!\n" );
  
     if ( Jet_BurstEndTime >= END_T )
       Aux_Error( ERROR_INFO, "Jet_BurstEndTime >= END_T !!\n" );
  
     if ( Flag_Burst4Vel && Jet_Burst4VelRatio <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_Burst4VelRatio <= Eps_double !!\n" );
  
     if ( Flag_BurstDens && Jet_BurstDensRatio <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_BurstDensRatio <= Eps_double !!\n" );
  
     if ( Flag_BurstTemp && Jet_BurstTempRatio <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_BurstTempRatio <= Eps_double !!\n" );
   }

   if ( IsothermalSlab_Center[0] == -1.0 )
        IsothermalSlab_Center[0] = 0.5*amr->BoxSize[0];

   if ( IsothermalSlab_Center[1] == -1.0 )
        IsothermalSlab_Center[1] = 0.5*amr->BoxSize[1];

   if ( IsothermalSlab_Center[2] == -1.0 )
        IsothermalSlab_Center[2] = 0.5*amr->BoxSize[2];

   if ( Jet_Ambient == 9 && OPT__INIT != 3 )
   {
      Aux_Error( ERROR_INFO, "OPT__INIT must be 3 !!\n" );
   }


// check UNIT_L is in reasonable range
   if ( ( UNIT_L <= 0.5*Const_kpc || 2.0*Const_kpc <= UNIT_L ) && OPT__UNIT )
      Aux_Error( ERROR_INFO, "UNIT_L=%e is far from %e !!\n", UNIT_L, Const_kpc );

   const double Const_Erg2eV = 6.2415e11;

// (1-2) convert to code unit
   Jet_SrcVel               *= Const_c     / UNIT_V;
   Jet_SrcTemp              *= Const_kB    / (ParticleMass*Const_c*Const_c);
   Jet_SrcDens              *= 1.0         / UNIT_D;

   Jet_Radius               *= Const_kpc   / UNIT_L;
   Jet_HalfHeight           *= Const_kpc   / UNIT_L;
   Jet_HalfOpeningAngle     *= M_PI        / 180.0;
   Jet_PrecessionAngle      *= M_PI        / 180.0;

   Jet_CenOffset[0]         *= Const_kpc   / UNIT_L;
   Jet_CenOffset[1]         *= Const_kpc   / UNIT_L;
   Jet_CenOffset[2]         *= Const_kpc   / UNIT_L;

   Jet_Duration             *= Const_Myr   / UNIT_T;

   if ( Jet_Ambient == 0  )
   {
     Amb_UniformDens        *= 1.0         / UNIT_D;
     Amb_UniformVel[0]      *= Const_c     / UNIT_V;
     Amb_UniformVel[1]      *= Const_c     / UNIT_V;
     Amb_UniformVel[2]      *= Const_c     / UNIT_V;
   }
   else if ( Jet_Ambient == 2 )
   {
     IsothermalSlab_VelocityDispersion  = Header[15];
     IsothermalSlab_PeakDens            = Header[16];
     //IsothermalSlab_Truncation          = Header[21];
     IsothermalSlab_Truncation          = 0.95*0.5*amr->BoxSize[2];
     gasDiskPeakDens                    = Header[19];

     IsothermalSlab_VelocityDispersion *= 1e5/UNIT_V; // km/s --> 1/c
     IsothermalSlab_PeakDens           *= 1.0         / UNIT_D;
     IsothermalSlab_Truncation         *= Const_kpc   / UNIT_L;
     IsothermalSlab_Center[0]          *= Const_kpc   / UNIT_L;
     IsothermalSlab_Center[1]          *= Const_kpc   / UNIT_L;
     IsothermalSlab_Center[2]          *= Const_kpc   / UNIT_L;
 
     distance_h                         = Header[29];
     v_halo                             = Header[30];

     distance_h                        *= Const_kpc   / UNIT_L;
     v_halo                            *= 1.0         / UNIT_V;
  

     gasDiskPeakDens                   /= UNIT_D;

     ambientTemperature                = Header[20]; 
     ambientTemperature               *= Const_kB/(ParticleMass*UNIT_V*UNIT_V);

     gasDiskTemperature                = Header[18];
     gasDiskTemperature               *= Const_kB/(ParticleMass*UNIT_V*UNIT_V);
   }
   

   Amb_UniformTemp          *= Const_kB    / (ParticleMass*Const_c*Const_c);
   Jet_AngularVelocity      *= 1.0;    // the unit of Jet_AngularVelocity is UNIT_T

   
   if ( Amb_FluSphereRadius > 0.0 )
   {
      Amb_FluSphereRadius   *= Const_kpc / UNIT_L;
   }


   if ( Jet_TimeDependentSrc )
   {
     Jet_BurstStartTime     *= 1e3*Const_yr / UNIT_T;
     Jet_BurstEndTime       *= 1e3*Const_yr / UNIT_T;
     Jet_Burst4VelRatio     *= Const_c      / UNIT_V;
     Jet_BurstDensRatio     *= 1.0          / UNIT_D;
   }


// (2) set the problem-specific derived parameters
   double SecAngle = 1.0 / cos(0.5*Jet_HalfOpeningAngle);
   double TanAngle = sin(0.5*Jet_HalfOpeningAngle) * SecAngle;

   Jet_MaxDis  = sqrt( SQR( Jet_Radius ) + SQR( Jet_HalfHeight * SecAngle ) + 2.0 * Jet_Radius * Jet_HalfHeight * TanAngle );

   for (int d=0; d<3; d++)    Jet_Center[d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[d];


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 0.5*BOX_SIZE * UNIT_L / (CharacteristicSpeed *UNIT_V) / UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      if ( CharacteristicSpeed == -1.0 ) Aux_Error( ERROR_INFO, "CharacteristicSpeed must be provided !!\n" );
      else                               END_T = End_T_Default;

      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

   if ( OUTPUT_DT < 0.0 )
   {
      OUTPUT_DT = END_T / 30.0;
      PRINT_WARNING( "OUTPUT_DT", OUTPUT_DT, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "=============================================================================\n" );
     Aux_Message( stdout, "  test problem ID          = %d\n",                TESTPROB_ID                                     );
     Aux_Message( stdout, "  Jet_Ambient              = %d\n",                Jet_Ambient                                     );
     Aux_Message( stdout, "  Jet_Fire                 = %d\n",                Jet_Fire                                        );
     Aux_Message( stdout, "  Jet_SmoothVel            = %d\n",                Jet_SmoothVel                                   );
     Aux_Message( stdout, "  Jet_Precession           = %d\n",                Jet_Precession                                  );
     Aux_Message( stdout, "  Jet_TimeDependentSrc     = %d\n",                Jet_TimeDependentSrc                            );
     Aux_Message( stdout, "  Jet_Duration             = %14.7e Myr \n",       Jet_Duration*UNIT_T/Const_Myr                   );
     Aux_Message( stdout, "  ParticleMass             = %14.7e g\n",          ParticleMass                                    );
     Aux_Message( stdout, "  Jet_SrcVel               = %14.7e c\n",          Jet_SrcVel                                      );
     Aux_Message( stdout, "  Jet_SrcDens              = %14.7e g/cm^3\n",     Jet_SrcDens*UNIT_D                              );
     Aux_Message( stdout, "  Jet_SrcTemp              = %14.7e kT/mc**2\n",   Jet_SrcTemp                                     );
     Aux_Message( stdout, "  Jet_NumDensSrc           = %14.7e per cc\n",     Jet_SrcDens*UNIT_D/ParticleMass                 );
     Aux_Message( stdout, "  Jet_CenOffset[x]         = %14.7e kpc\n",        Jet_CenOffset [0]*UNIT_L/Const_kpc              );
     Aux_Message( stdout, "  Jet_CenOffset[y]         = %14.7e kpc\n",        Jet_CenOffset [1]*UNIT_L/Const_kpc              );
     Aux_Message( stdout, "  Jet_CenOffset[z]         = %14.7e kpc\n",        Jet_CenOffset [2]*UNIT_L/Const_kpc              );
     Aux_Message( stdout, "  Jet_AngularVelocity      = %14.7e degree/kyr\n", Jet_AngularVelocity                             );
     Aux_Message( stdout, "  Jet_PrecessionAngle      = %14.7e degree\n",     Jet_PrecessionAngle*180.0/M_PI                  );
     Aux_Message( stdout, "  Jet_HalfOpeningAngle     = %14.7e degree\n",     Jet_HalfOpeningAngle*180.0/M_PI                 );
     Aux_Message( stdout, "  Jet_Radius               = %14.7e kpc\n",        Jet_Radius*UNIT_L/Const_kpc                     );
     Aux_Message( stdout, "  Jet_HalfHeight           = %14.7e kpc\n",        Jet_HalfHeight*UNIT_L/Const_kpc                 );
     Aux_Message( stdout, "  Jet_MaxDis               = %14.7e kpc\n",        Jet_MaxDis*UNIT_L/Const_kpc                     );
   }

   if ( Jet_Ambient == 0 && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Amb_UniformDens          = %14.7e g/cm^3\n",     Amb_UniformDens*UNIT_D                          );
     Aux_Message( stdout, "  Amb_UniformTemp          = %14.7e kT/mc**2\n",   Amb_UniformTemp                                 );
     Aux_Message( stdout, "  Amb_UniformVel[x]        = %14.7e c\n",          Amb_UniformVel[0]                               );
     Aux_Message( stdout, "  Amb_UniformVel[y]        = %14.7e c\n",          Amb_UniformVel[1]                               );
     Aux_Message( stdout, "  Amb_UniformVel[z]        = %14.7e c\n",          Amb_UniformVel[2]                               );
     Aux_Message( stdout, "  Jet_UniformNumDens       = %14.7e per cc\n",     Amb_UniformDens*UNIT_D/ParticleMass             );
   }
   else if ( Jet_Ambient == 2 && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  IsothermalSlab_Center[0]          = %14.7e kpc\n",      IsothermalSlab_Center[0]*UNIT_L/Const_kpc );
     Aux_Message( stdout, "  IsothermalSlab_Center[1]          = %14.7e kpc\n",      IsothermalSlab_Center[1]*UNIT_L/Const_kpc );
     Aux_Message( stdout, "  IsothermalSlab_Center[2]          = %14.7e kpc\n",      IsothermalSlab_Center[2]*UNIT_L/Const_kpc );
   }

   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  CharacteristicSpeed      = %14.7e c\n",          CharacteristicSpeed / UNIT_V                    );
   }


   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_PrecessionAxis[x]    = %14.7e\n",            Jet_PrecessionAxis[0]                           );
     Aux_Message( stdout, "  Jet_PrecessionAxis[y]    = %14.7e\n",            Jet_PrecessionAxis[1]                           );
     Aux_Message( stdout, "  Jet_PrecessionAxis[z]    = %14.7e\n",            Jet_PrecessionAxis[2]                           );
   }

   if ( Amb_FluSphereRadius > 0.0 && MPI_Rank == 0 ) 
   {
     Aux_Message( stdout, "  Amb_FluSphereRadius      = %14.7e kpc\n",        Amb_FluSphereRadius*UNIT_L/Const_kpc           );

   }

   if ( Jet_TimeDependentSrc && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_BurstStartTime       = %14.7e kyr \n",       Jet_BurstStartTime*UNIT_T/(1e3*Const_yr)        );
     Aux_Message( stdout, "  Jet_BurstEndTime         = %14.7e kyr \n",       Jet_BurstEndTime*UNIT_T/(1e3*Const_yr)          );
     Aux_Message( stdout, "  Jet_Burst4VelRatio       = %14.7e c \n",         Jet_Burst4VelRatio                              );
     Aux_Message( stdout, "  Jet_BurstDensRatio       = %14.7e g/cm^3\n",     Jet_BurstDensRatio*UNIT_D                       );
     Aux_Message( stdout, "  Jet_BurstTempRatio       = %14.7e\n",            Jet_BurstTempRatio                              );
     Aux_Message( stdout, "  Flag_Burst4Vel           = %d\n",                Flag_Burst4Vel                                  );
     Aux_Message( stdout, "  Flag_BurstDens           = %d\n",                Flag_BurstDens                                  );
     Aux_Message( stdout, "  Flag_BurstTemp           = %d\n",                Flag_BurstTemp                                  );
   }

   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "=============================================================================\n"                   );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n"                                    );

} // FUNCTION : SetParameter


void ReadBinFile( char *FileName, real **buffer)
{
  FILE *pFile;
  long lSize;
  size_t result;

  pFile = fopen ( FileName , "rb" );
  if (pFile==NULL) {fputs ("File error\n",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  *buffer = (real*) calloc (lSize,sizeof(double));
  if (*buffer == NULL) {fputs ("Memory error\n",stderr); exit (2);}

  // copy the file into the *buffer:
  result = fread (*buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error\n",stderr); exit (3);}

  fclose (pFile);
}


void SetArray()
{
// Reading table for interpolations in SetGridIC()
   char TableFileName[] = "UM_IC";
   ReadBinFile(TableFileName, &buffer);

   int headerSize = (int)buffer[0];

   Header = (real*)malloc((size_t)headerSize*sizeof(real));

   memcpy(Header, buffer, (size_t)headerSize*sizeof(real));

   ParticleMass = Header[8]*Header[9];

   int Nx = (int)Header[22];
   int Ny = (int)Header[23];
   int Nz = (int)Header[24];

   if (5*Nx*Ny*Nz > INT_MAX) {printf("integer overflow!!\n"); exit(0);}

   real Lx = Header[12];
   real Ly = Header[13];
   real Lz = Header[14];
   int numGhost = (int)Header[25];

   real dx = Lx/(real)Nx;
   real dy = Lx/(real)Nx;
   real dz = Lx/(real)Nx;

   int NX = Nx+2*numGhost;
   int NY = Ny+2*numGhost;
   int NZ = Nz+2*numGhost;

   Rhoo = (real***)calloc_3d_array((size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real));
   VelX = (real***)calloc_3d_array((size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real));
   VelY = (real***)calloc_3d_array((size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real));
   VelZ = (real***)calloc_3d_array((size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real));
   Pres = (real***)calloc_3d_array((size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real));

   X = (real*)calloc((size_t)NX,sizeof(real));
   Y = (real*)calloc((size_t)NY,sizeof(real));
   Z = (real*)calloc((size_t)NZ,sizeof(real));


   //for (int i=-numGhost;i<NX-numGhost;i++) X[i+numGhost] = (0.5+(real)i)*dx;
   //for (int i=-numGhost;i<NY-numGhost;i++) Y[i+numGhost] = (0.5+(real)i)*dy;
   //for (int i=-numGhost;i<NZ-numGhost;i++) Z[i+numGhost] = (0.5+(real)i)*dz;

   real *Ptr;

   Ptr = buffer + headerSize;


   for (int c=0;c<5*NX*NY*NZ;c++){
     int i, j, k, cc;

     cc = c%(NX*NY*NZ);
     i = (cc - cc%(NY*NZ)) / (NY*NZ);
     j = ((cc - cc%NZ) / NZ) % NY;
     k = cc%NZ;

     if (          0 <= c && c <   NX*NY*NZ ) Rhoo[i][j][k] = Ptr[c];
     if (   NX*NY*NZ <= c && c < 2*NX*NY*NZ ) VelX[i][j][k] = Ptr[c];
     if ( 2*NX*NY*NZ <= c && c < 3*NX*NY*NZ ) VelY[i][j][k] = Ptr[c];
     if ( 3*NX*NY*NZ <= c && c < 4*NX*NY*NZ ) VelZ[i][j][k] = Ptr[c];
     if ( 4*NX*NY*NZ <= c && c < 5*NX*NY*NZ ) Pres[i][j][k] = Ptr[c];

   }
  
   Ptr += 5*NX*NY*NZ;
   for (int c=0;c<NX;c++) X[c] = Ptr[c];

   Ptr += NX;
   for (int c=0;c<NY;c++) Y[c] = Ptr[c];

   Ptr += NY;
   for (int c=0;c<NZ;c++) Z[c] = Ptr[c];

   //bool Pass = true;

   //for (int i=0;i<NX;i++){
   //for (int j=0;j<NY;j++){
   //for (int k=0;k<NZ;k++){
   //
   //Pass &= Rhoo[i][j][k] == (real)16.;
   //Pass &= VelX[i][j][k] == (real)32.;
   //Pass &= VelY[i][j][k] == (real)64.;
   //Pass &= VelZ[i][j][k] == (real)128.;
   //Pass &= Pres[i][j][k] == (real)256.;

   //}}}
   //if (Pass == false){ printf("fail!!!!!!!!!!!!!!!!!!!\n"); exit(0); }
   //if (Pass == true) { printf("pass!!!!!!!!!!!!!!!!!!!\n"); exit(0); }
}


void Interpolation_UM_IC( real x, real y, real z, real *Pri )
{
  real xyz[3] = {x, y, z};
  int Nx = (int)Header[22];
  int Ny = (int)Header[23];
  int Nz = (int)Header[24];

  real dx = Header[26];
  real dy = Header[27];
  real dz = Header[28];

  real dxyz[3] = {dx, dy, dz};
 
  int numGhost = (int)Header[25];

  int NX = Nx+2*numGhost;
  int NY = Ny+2*numGhost;
  int NZ = Nz+2*numGhost;

  int Idx = Mis_BinarySearch_Real(X, 0, NX-1, x);
  int Jdx = Mis_BinarySearch_Real(Y, 0, NY-1, y);
  int Kdx = Mis_BinarySearch_Real(Z, 0, NZ-1, z);

  if (Idx<0 || Idx > NX-2){ printf("x=%e is out of range! X[0]=%e, X[%d]=%e\n", x, X[0], NX-1, X[NX-1]); exit(0); }
  if (Jdx<0 || Jdx > NY-2){ printf("y=%e is out of range! Y[0]=%e, Y[%d]=%e\n", y, Y[0], NY-1, Y[NY-1]); exit(0); }
  if (Kdx<0 || Kdx > NZ-2){ printf("z=%e is out of range! Z[0]=%e, Z[%d]=%e\n", z, Z[0], NZ-1, Z[NZ-1]); exit(0); }


  real Vertex000[5] = {Rhoo[Idx  ][Jdx  ][Kdx  ], VelX[Idx  ][Jdx  ][Kdx  ], VelY[Idx  ][Jdx  ][Kdx  ], VelZ[Idx  ][Jdx  ][Kdx  ], Pres[Idx  ][Jdx  ][Kdx  ]};
  real Vertex001[5] = {Rhoo[Idx  ][Jdx  ][Kdx+1], VelX[Idx  ][Jdx  ][Kdx+1], VelY[Idx  ][Jdx  ][Kdx+1], VelZ[Idx  ][Jdx  ][Kdx+1], Pres[Idx  ][Jdx  ][Kdx+1]};
  real Vertex010[5] = {Rhoo[Idx  ][Jdx+1][Kdx  ], VelX[Idx  ][Jdx+1][Kdx  ], VelY[Idx  ][Jdx+1][Kdx  ], VelZ[Idx  ][Jdx+1][Kdx  ], Pres[Idx  ][Jdx+1][Kdx  ]};
  real Vertex100[5] = {Rhoo[Idx+1][Jdx  ][Kdx  ], VelX[Idx+1][Jdx  ][Kdx  ], VelY[Idx+1][Jdx  ][Kdx  ], VelZ[Idx+1][Jdx  ][Kdx  ], Pres[Idx+1][Jdx  ][Kdx  ]};
  real Vertex011[5] = {Rhoo[Idx  ][Jdx+1][Kdx+1], VelX[Idx  ][Jdx+1][Kdx+1], VelY[Idx  ][Jdx+1][Kdx+1], VelZ[Idx  ][Jdx+1][Kdx+1], Pres[Idx  ][Jdx+1][Kdx+1]};
  real Vertex101[5] = {Rhoo[Idx+1][Jdx  ][Kdx+1], VelX[Idx+1][Jdx  ][Kdx+1], VelY[Idx+1][Jdx  ][Kdx+1], VelZ[Idx+1][Jdx  ][Kdx+1], Pres[Idx+1][Jdx  ][Kdx+1]};
  real Vertex110[5] = {Rhoo[Idx+1][Jdx+1][Kdx  ], VelX[Idx+1][Jdx+1][Kdx  ], VelY[Idx+1][Jdx+1][Kdx  ], VelZ[Idx+1][Jdx+1][Kdx  ], Pres[Idx+1][Jdx+1][Kdx  ]};
  real Vertex111[5] = {Rhoo[Idx+1][Jdx+1][Kdx+1], VelX[Idx+1][Jdx+1][Kdx+1], VelY[Idx+1][Jdx+1][Kdx+1], VelZ[Idx+1][Jdx+1][Kdx+1], Pres[Idx+1][Jdx+1][Kdx+1]};

  bool Unphy = false;

  Unphy |= SRHD_CheckUnphysical( NULL, Vertex000, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  );
  Unphy |= SRHD_CheckUnphysical( NULL, Vertex001, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  );
  Unphy |= SRHD_CheckUnphysical( NULL, Vertex010, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  );
  Unphy |= SRHD_CheckUnphysical( NULL, Vertex100, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  );
  Unphy |= SRHD_CheckUnphysical( NULL, Vertex011, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  );
  Unphy |= SRHD_CheckUnphysical( NULL, Vertex110, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  );
  Unphy |= SRHD_CheckUnphysical( NULL, Vertex101, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  );
  Unphy |= SRHD_CheckUnphysical( NULL, Vertex111, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  );

  if (Unphy){
    printf("Idx=%d, Jdx=%d, Kdx=%d\n", Idx, Jdx, Kdx);
    printf("x=%e, y=%e, z=%e\n", x, y, z);
    exit(0);
  }

  real xyz000[3] = {X[Idx], Y[Jdx], Z[Kdx]};

  for(int v=0;v<5;v++){
    real FieldAtVertices[8] = {Vertex000[v], Vertex001[v], Vertex010[v], Vertex100[v], Vertex011[v], Vertex101[v], Vertex110[v], Vertex111[v]};

    Pri[v] = TrilinearInterpolation(FieldAtVertices, xyz000, dxyz, xyz);
  }

  //free_3d_array((void***)Rhoo);
  //free_3d_array((void***)VelX);
  //free_3d_array((void***)VelY);
  //free_3d_array((void***)VelZ);
  //free_3d_array((void***)Pres);
  //free(X);
  //free(Y);
  //free(Z);
  //free(buffer);
}

#ifdef GRAVITY
real IsothermalSlab_Pot(real z)
{
  real Pot, Log;

  // 1. isothermal slab
  Pot  = 2.0*M_PI*NEWTON_G*IsothermalSlab_PeakDens;
  Pot /= SQR(IsothermalSlab_VelocityDispersion);
  Pot  = log(cosh(z*sqrt(Pot)));
  Pot *= 2.0*SQR(IsothermalSlab_VelocityDispersion);
  
  // 2. log potential
  Log  = SQR(v_halo) * log(z*z + SQR(distance_h));
  Log -= SQR(v_halo) * log(SQR(distance_h));

  return Pot + Log;
}
#endif

//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{
// variables for jet
   real Pri[NCOMP_FLUID];
   real xc = x - IsothermalSlab_Center[0];
   real yc = y - IsothermalSlab_Center[1];
   real zc = z - IsothermalSlab_Center[2];
   interfaceHeight = Header[17];
   interfaceHeight *= Const_kpc/UNIT_L;

   if ( Jet_Ambient == 0 ) // uniform ambient
   {
      Pri[0] = (real)Amb_UniformDens;
      Pri[1] = (real)Amb_UniformVel[0];
      Pri[2] = (real)Amb_UniformVel[1];
      Pri[3] = (real)Amb_UniformVel[2];
      Pri[4] = (real)Amb_UniformTemp * Amb_UniformDens;
   }
   else if ( Jet_Ambient == 2 )
   {
#    ifdef GRAVITY
     if (fabs(zc) < interfaceHeight){
       Interpolation_UM_IC( xc, yc, zc, Pri);

       if ( SRHD_CheckUnphysical( NULL, Pri,
                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,                              
                                  EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  ) ) exit(0);
     }
     else{

       real Dens_gDisk_ambient, PotAtZ0, ambientDens;

       PotAtZ0 = IsothermalSlab_Pot(interfaceHeight);

       Dens_gDisk_ambient  = ambientTemperature / gasDiskTemperature;
       Dens_gDisk_ambient *= exp( PotAtZ0*(ambientTemperature-gasDiskTemperature)/(ambientTemperature*gasDiskTemperature) );

       if (Dens_gDisk_ambient > HUGE_NUMBER || Dens_gDisk_ambient < -HUGE_NUMBER){
         printf("Dens_gDisk_ambient=%e! %s: %d\n", Dens_gDisk_ambient, __FUNCTION__, __LINE__);
         exit(0);
       }

       real ambientPeakDens  = gasDiskPeakDens / Dens_gDisk_ambient; 
      
       if (fabs(zc) > IsothermalSlab_Truncation)
           ambientDens  = -IsothermalSlab_Pot(IsothermalSlab_Truncation)/ambientTemperature;
       else 
           ambientDens  = -IsothermalSlab_Pot(zc)/ambientTemperature;

       ambientDens  = exp(ambientDens);
       ambientDens *= ambientPeakDens;
   
       Pri[0] = ambientDens;
       Pri[1] = 0.0;
       Pri[2] = 0.0;
       Pri[3] = 0.0;
       Pri[4] = ambientDens*ambientTemperature;

       if ( SRHD_CheckUnphysical( NULL, Pri,
                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,                              
                                  EoS_AuxArray_Int, h_EoS_Table,  __FUNCTION__, __LINE__, true  ) ) exit(0);
     }
#  endif
   }


   Hydro_Pri2Con( Pri, fluid, NULL_BOOL, NULL_INT, NULL, EoS_DensPres2Eint_CPUPtr,
                  EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

} // FUNCTION : SetGridIC






//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Jets
// Description :  Function to reset the fluid field
//
// Note        :  1. Invoked by "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()" using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                2. This function will be invoked when constructing the initial condition
//                    (by calling "Model_Init_ByFunction_AssignData()") and after each update
//                    (by calling "Flu_ResetByUser_API()")
//                3. Input "fluid" array stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()"
//                5. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//
// =======================================================================================
/*        G       A       C              */ 
/*          ____________                 */
/*          \     |     /                */
/*           \   E|    /        z        */
/*            \  /|\  /         ^        */
/*             \/_|_\/          |        */
/*             /\O| /\B         |        */
/*            /  \|/  \                  */
/*           /   D|    \                 */
/*          /_____|_____\                */
/*                F                      */
// =======================================================================================
//
bool Flu_ResetByUser_Jets( real fluid[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] )
{
  if ( Jet_Fire == 0 ) return false;

  if ( Jet_Duration < Time ) return false;

  double xp[3], rp[3];
  double Prim[5], Cons[5], Vel[3];
  real PriReal[5];
  double PrecessionAxis_Spherical[3], Omega_t;
  bool InsideUpperCone, InsideLowerCone;
  double Jet_SrcVelSmooth;

  Omega_t = Jet_AngularVelocity * Time * M_PI / 180.0;



  // shift the coordinate origin to the source center (the point O)
  xp[0] = x - Jet_Center[0];
  xp[1] = y - Jet_Center[1];
  xp[2] = z - Jet_Center[2];

  if ( Jet_PrecessionAxis[0] != 0.0 || Jet_PrecessionAxis[1] != 0.0 ||  Jet_PrecessionAxis[2] == 0.0 )
  {
    // get theta, phi for the first rotation
    Mis_Cartesian2Spherical( Jet_PrecessionAxis, PrecessionAxis_Spherical );

    // rotate coordinate to align z-axis with fixed precession axis
	CartesianRotate(xp, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], false);
  }

  // rotate coordinate to align z-axis with rotating symmetric axis
  CartesianRotate(xp, Jet_PrecessionAngle, Omega_t, false);

  // determine whether or not the point is inside of source
  InsideUpperCone  = SQR(xp[0]) + SQR(xp[1]) <= SQR( +tan(Jet_HalfOpeningAngle)*xp[2] + Jet_Radius );
  InsideUpperCone &= 0.0 <= xp[2] && xp[2] <= Jet_HalfHeight;

  InsideLowerCone  = SQR(xp[0]) + SQR(xp[1]) <= SQR( -tan(Jet_HalfOpeningAngle)*xp[2] + Jet_Radius );
  InsideLowerCone &= -Jet_HalfHeight <= xp[2] && xp[2] <= 0.0;


  if ( Jet_HalfOpeningAngle != 0.0 )
  {
   InsideUpperCone &= SQR(xp[0]) + SQR(xp[1]) + SQR(xp[2] + Jet_Radius/tan(Jet_HalfOpeningAngle)) 
                   <= SQR(Jet_HalfHeight+Jet_Radius/tan(Jet_HalfOpeningAngle));

   InsideUpperCone &= 0.0 <= xp[2]; 

   InsideLowerCone &= SQR(xp[0]) + SQR(xp[1]) + SQR(xp[2] - Jet_Radius/tan(Jet_HalfOpeningAngle)) 
                   <= SQR(Jet_HalfHeight+Jet_Radius/tan(Jet_HalfOpeningAngle));

   InsideLowerCone &= xp[2] <= 0.0;
  }
  else
  {
    InsideUpperCone &= 0.0 <= xp[2] && xp[2] <= Jet_HalfHeight;
    InsideLowerCone &= -Jet_HalfHeight <= xp[2] && xp[2] <= 0.0;
  }



  // set fluid variable inside source
  if ( ( InsideUpperCone && ( Jet_Fire == 1 || Jet_Fire == 3 ) ) 
	|| ( InsideLowerCone && ( Jet_Fire == 2 || Jet_Fire == 3 ) ) )
  {
    if ( Jet_HalfOpeningAngle == 0.0 )
  	{
  	  Vel[0] = 0.0;
  	  Vel[1] = 0.0;

  	  if ( InsideUpperCone == true ) Vel[2] = +Jet_SrcVel;
      else                           Vel[2] = -Jet_SrcVel;

  	  CartesianRotate(Vel, Jet_PrecessionAngle, Omega_t, true);

      if ( Jet_PrecessionAxis[0] != 0.0 || Jet_PrecessionAxis[1] != 0.0 ||  Jet_PrecessionAxis[2] == 0.0 )
           CartesianRotate(Vel, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], true);

      Prim[0] = Jet_SrcDens;
  	  Prim[1] = Vel[0];
  	  Prim[2] = Vel[1];
  	  Prim[3] = Vel[2];
      Prim[4] = Jet_SrcTemp*Jet_SrcDens;
  	}
  	else
  	{
      // shift origin to the point D/E
	  if ( InsideUpperCone == true )   xp[2] += Jet_Radius/tan(Jet_HalfOpeningAngle);
	  else                             xp[2] -= Jet_Radius/tan(Jet_HalfOpeningAngle);

      CartesianRotate(xp, Jet_PrecessionAngle, Omega_t, true);

      Mis_Cartesian2Spherical(xp, rp);

      if ( InsideLowerCone == true ) rp[1] -= M_PI;

      // smooth velocity on cross section
      if ( Jet_SmoothVel ) Jet_SrcVelSmooth = Jet_SrcVel*SQR(cos( 0.5 * M_PI * rp[1] / Jet_HalfOpeningAngle ));
      else                 Jet_SrcVelSmooth = Jet_SrcVel;


      if ( Jet_PrecessionAxis[0] != 0.0 || Jet_PrecessionAxis[1] != 0.0 ||  Jet_PrecessionAxis[2] == 0.0 )
         CartesianRotate(xp, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], true);

      Mis_Cartesian2Spherical(xp, rp);

      Prim[0] = Jet_SrcDens;
      Prim[1] = Jet_SrcVelSmooth*sin(rp[1])*cos(rp[2]);
      Prim[2] = Jet_SrcVelSmooth*sin(rp[1])*sin(rp[2]);
      Prim[3] = Jet_SrcVelSmooth*cos(rp[1]);
      Prim[4] = Jet_SrcTemp*Jet_SrcDens;
  	}

    PriReal[0] = (real)Prim[0];
    PriReal[1] = (real)Prim[1];
    PriReal[2] = (real)Prim[2];
    PriReal[3] = (real)Prim[3];
    PriReal[4] = (real)Prim[4];

    Hydro_Pri2Con( PriReal, fluid, NULL_BOOL, NULL_INT, NULL, EoS_DensPres2Eint_CPUPtr,
                   EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                   EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

	return true;
  }



  return false;


} // FUNCTION : Flu_ResetByUser_Jets


// (true/false): if the target cell (is/is not) within the region to be refined
static bool Flag_Region( const int i, const int j, const int k, const int lv, const int PID )
{

    const double dh     = amr->dh[lv];                                                         // grid size
    const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,  // x,y,z position
                            amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                            amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   
    bool Flag = false;  
   
    const double Center[3]      = { 0.5*amr->BoxSize[0], 
                                    0.5*amr->BoxSize[1], 
                                    0.5*amr->BoxSize[2] };
   
    //const double dR[3]          = { Pos[0]-Center[0]-Jet_CenOffset[0], 
    //                                Pos[1]-Center[1]-Jet_CenOffset[1], 
    //                                Pos[2]-Center[2]-Jet_CenOffset[2] };
   
    //const double R              = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );
   
    //const double ShellThickness = 16*amr->dh[3];
    
   
   
    if ( fabs(Pos[2]-Center[2]) <= interfaceHeight )  return true;
    else if ( lv >= 4)                                return false;
    else                                              return true;


} // FUNCTION : Flag_Region



bool Flag_User( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{
  const double dh     = amr->dh[lv];                                                  // grid size
  const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                          amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                          amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

  const double Center[3]      = { Jet_Center[0], 
                                  Jet_Center[1], 
                                  Jet_Center[2] };

  const double dR[3]          = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
  const double R              = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

  
  bool Flag;


  Flag  = fabs(Pos[2] - Center[2]) < dh*1.8;

  bool Disk, Bubble, Src;

  //Disk     = fabs(Pos[2]-Center[2]) <= interfaceHeight;
  //Disk    &= R > Jet_Radius;
  //Disk    &= lv <= Disk_MAX_LEVEL;

  //Bubble   = fabs(Pos[2]-Center[2]) > interfaceHeight;
  //Bubble  &= lv <= Bubb_MAX_LEVEL;
  
  //Src      = fabs(Pos[2]-Center[2]) <= interfaceHeight;
  //Src     &= R <= Jet_Radius;

  Src      = R <= dh*1.8;

  //return Disk || Bubble || Src;
  return Src;
} // FUNCTION : Flag_User

void CartesianRotate( double x[], double theta, double phi, bool inverse )
{
  double xp[3];

  if ( inverse )
  {
     xp[0] = -            sin(phi)*x[0] - cos(theta)*cos(phi)*x[1] + sin(theta)*cos(phi)*x[2];
	 xp[1] = +            cos(phi)*x[0] - cos(theta)*sin(phi)*x[1] + sin(theta)*sin(phi)*x[2];
	 xp[2] =                            + sin(theta)*         x[1] + cos(theta)*         x[2];
  }
  else
  {
     xp[0] = -            sin(phi)*x[0] +            cos(phi)*x[1];
	 xp[1] = - cos(theta)*cos(phi)*x[0] - cos(theta)*sin(phi)*x[1] + sin(theta)*         x[2];
	 xp[2] = + sin(theta)*cos(phi)*x[0] + sin(theta)*sin(phi)*x[1] + cos(theta)*         x[2];
  }

  for (int i=0;i<3;i++) x[i] = xp[i];
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Jets
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Jets()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// set the problem-specific runtime parameters
   SetParameter();


// get enclosed mass
   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = Flag_User;
   Flag_Region_Ptr          = Flag_Region;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser_Jets;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;
#  ifdef GRAVITY
   Init_ExtPot_Ptr          = Init_ExtPot_IsothermalSlab;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Jets

#endif
