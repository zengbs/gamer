#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#endif

#ifdef GRAVITY




// =================================
// I. Set auxiliary arrays
// =================================

#ifndef __CUDACC__

extern double MilkyWay_Halo_v;
extern double MilkyWay_Halo_d;
extern double MilkyWay_Disk_M;
extern double MilkyWay_Disk_a;
extern double MilkyWay_Disk_b;
extern double MilkyWay_Bulge_M;
extern double MilkyWay_Bulge_d;
extern double MilkyWay_Temperature;
extern double MilkyWay_Center[3];
extern int    MilkyWay_Trun;
extern double MilkyWay_TrunRhoRatio;

//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtPotAuxArray_MilkyWay
// Description :  Set the auxiliary arrays ExtPot_AuxArray_Flt/Int[] used by ExtPot_MilkyWay()
//
// Note        :  1. Invoked by Init_ExtPot_MilkyWay()
//                2. AuxArray_Flt/Int[] have the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
void SetExtPotAuxArray_MilkyWay( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[ 0] = MilkyWay_Center[0];          // x coordinate of the MilkyWay center
   AuxArray_Flt[ 1] = MilkyWay_Center[1];          // y ...
   AuxArray_Flt[ 2] = MilkyWay_Center[2];          // z ...
   AuxArray_Flt[ 3] = SQR(MilkyWay_Halo_v);        // 
   AuxArray_Flt[ 4] = SQR(MilkyWay_Halo_d);        // 
   AuxArray_Flt[ 5] = MilkyWay_Disk_a;             // 
   AuxArray_Flt[ 6] = SQR(MilkyWay_Disk_b);        // 
   AuxArray_Flt[ 7] = MilkyWay_Bulge_d;            // 
   AuxArray_Flt[ 8] = NEWTON_G*MilkyWay_Disk_M;    // -G*M_disk
   AuxArray_Flt[ 9] = NEWTON_G*MilkyWay_Bulge_M;   // -G*M_bulge
   AuxArray_Flt[10] = MilkyWay_TrunRhoRatio;
   AuxArray_Flt[11] = MilkyWay_Temperature;

   AuxArray_Int[ 0] = MilkyWay_Trun;

} // FUNCTION : SetExtPotAuxArray_MilkyWay
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external potential
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtPot_MilkyWay
// Description :  Calculate the external potential at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary arrays UserArray_Flt/Int[] are set by SetExtPotAuxArray_MilkyWay()
//
// Parameter   :  x/y/z             : Target spatial coordinates
//                Time              : Target physical time
//                UserArray_Flt/Int : User-provided floating-point/integer auxiliary arrays
//                Usage             : Different usages of external potential when computing total potential on level Lv
//                                    --> EXT_POT_USAGE_ADD     : add external potential on Lv
//                                        EXT_POT_USAGE_SUB     : subtract external potential for preparing self-gravity potential on Lv-1
//                                        EXT_POT_USAGE_SUB_TINT: like SUB but for temporal interpolation
//                                    --> This parameter is useless in most cases
//                PotTable          : 3D potential table used by EXT_POT_TABLE
//
// Return      :  External potential at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real ExtPot_MilkyWay( const double x, const double y, const double z, const double Time,
                             const double UserArray_Flt[], const int UserArray_Int[],
                             const ExtPotUsage_t Usage, const real PotTable[] )
{

// halo potential
   const double cx                    =       UserArray_Flt[ 0];   // x coordinate of the MilkyWay center
   const double cy                    =       UserArray_Flt[ 1];   // y ...
   const double cz                    =       UserArray_Flt[ 2];   // z ...
   const real   Halo_v2               = (real)UserArray_Flt[ 3];   // 
   const real   Halo_d2               = (real)UserArray_Flt[ 4];   // 
   const real   Disk_a                = (real)UserArray_Flt[ 5];   // 
   const real   Disk_b2               = (real)UserArray_Flt[ 6];   // 
   const real   Bulge_d               = (real)UserArray_Flt[ 7];   // 
   const real   GDiskM                = (real)UserArray_Flt[ 8];   // 
   const real   GBulgeM               = (real)UserArray_Flt[ 9];   // 
   const real   MilkyWay_TrunRhoRatio = (real)UserArray_Flt[10];   // 
   const real   MilkyWay_Temperature  = (real)UserArray_Flt[11];   // 

   const int   MilkyWay_Trun          =       UserArray_Int[ 0];   // 

   const real   dx  = (real)(x - cx);
   const real   dy  = (real)(y - cy);
   const real   dz  = (real)(z - cz);

   const real   r2  = SQR(dx) + SQR(dy) + SQR(dz);
   const real   r   = SQRT(r2);
   const real   dx2 = dx*dx;
   const real   dy2 = dy*dy;
   const real   dz2 = dz*dz;

   real   HaloPot  = Halo_v2*log( r2 + Halo_d2 );
   real   DiskPot  = -GDiskM;
          DiskPot /= SQRT( dx2 + dy2 + SQR( Disk_a + SQRT( dz2 + Disk_b2 ) ) );
   real   BulgePot = -GBulgeM/( r + Bulge_d );


   real TotPot        = HaloPot + DiskPot + BulgePot;
   real TotPotCenter  = Halo_v2*log( Halo_d2 );
        TotPotCenter -= GDiskM / SQRT( SQR( Disk_a + SQRT( Disk_b2 ) ) );
        TotPotCenter -= GBulgeM/ Bulge_d;

   if ( MilkyWay_Trun == 1 )
   {
      if ( TotPot - TotPotCenter > MilkyWay_Temperature*log(MilkyWay_TrunRhoRatio) )
           TotPot = TotPotCenter + MilkyWay_Temperature*log(MilkyWay_TrunRhoRatio);
   }

   return TotPot;

} // FUNCTION : ExtPot_MilkyWay



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtPot_t ExtPot_Ptr = ExtPot_MilkyWay;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtPot_MilkyWay
// Description :  Return the function pointers of the CPU/GPU external potential routines
//
// Note        :  1. Invoked by Init_ExtPot_MilkyWay()
//
// Parameter   :  CPU/GPUExtPot_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtPot_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtPot_MilkyWay( ExtPot_t &GPUExtPot_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtPot_Ptr, ExtPot_Ptr, sizeof(ExtPot_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtPot_MilkyWay( ExtPot_t &CPUExtPot_Ptr )
{
   CPUExtPot_Ptr = ExtPot_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtPotAuxArray_MilkyWay( double [], int [] );
void SetCPUExtPot_MilkyWay( ExtPot_t & );
#ifdef GPU
void SetGPUExtPot_MilkyWay( ExtPot_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtPot_MilkyWay
// Description :  Initialize external potential
//
// Note        :  1. Set auxiliary arrays by invoking SetExtPotAuxArray_*()
//                   --> They will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external potential major routines by invoking SetCPU/GPUExtPot_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtPot_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtPot_MilkyWay()
{

   SetExtPotAuxArray_MilkyWay( ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int );
   SetCPUExtPot_MilkyWay( CPUExtPot_Ptr );
#  ifdef GPU
   SetGPUExtPot_MilkyWay( GPUExtPot_Ptr );
#  endif

} // FUNCTION : Init_ExtPot_MilkyWay

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
