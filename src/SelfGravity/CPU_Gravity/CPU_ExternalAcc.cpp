#ifndef __EXTERNALACC__
#define __EXTERNALACC__



#include "CUPOT.h"

#ifdef GRAVITY



// soften length implementation
#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT




//-----------------------------------------------------------------------------------------
// Function    :  ExternalAcc
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                3. Axiliary array UserArray[] is set by "Init_ExternalAcc_Ptr", which
//                   points to Init_ExternalAcc() by default but may be overwritten by various
//                   test problem initializers
//                4. By default we assume
//                     UserArray[0] = x coordinate of the external acceleration center
//                     UserArray[1] = y ...
//                     UserArray[2] = z ..
//                     UserArray[3] = gravitational_constant*point_source_mass
//                     UserArray[4] = soften_length (<=0.0 --> disable)
//                   --> But one can easily modify this file to change the default behavior
//                5. Two different soften length implementations are supported
//                   --> SOFTEN_PLUMMER & SOFTEN_RUFFERT
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Current physical time
//                UserArray : User-provided auxiliary array (set by "Init_ExternalAcc_Ptr")
//
// Return      :  Acc
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void ExternalAcc( real Acc[], const double x, const double y, const double z, 
				  const double Time, ExtAcc_AuxStruct_t ExtAcc_AuxStruct )
{

   const real GM       = (real)ExtAcc_AuxStruct.GM;
   const real eps      = (real)ExtAcc_AuxStruct.Eps;
   const real dx       = (real)(x - ExtAcc_AuxStruct.Center[0]);
   const real dy       = (real)(y - ExtAcc_AuxStruct.Center[1]);
   const real dz       = (real)(z - ExtAcc_AuxStruct.Center[2]);
   const real r        = SQRT( dx*dx + dy*dy + dz*dz );

// Plummer
#  if   ( defined SOFTEN_PLUMMER )
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps), (real)-1.5 );

// Ruffert 1994
#  elif ( defined SOFTEN_RUFFERT )
   const real tmp = EXP( -SQR(r)/SQR(eps) );
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps)*tmp, (real)-1.5 )*( (real)1.0 - tmp );

#  else
   const real _r3 = (real)1.0/CUBE(r);
#  endif

   Acc[0] = -GM*_r3*dx;
   Acc[1] = -GM*_r3*dy;
   Acc[2] = -GM*_r3*dz;

   real *Table_Acc, *Table_R;
   int ExtAcc_NBin;

   Table_Acc    =ExtAcc_AuxStruct.ExtAcc_Table_Acc;
   Table_R     = ExtAcc_AuxStruct.ExtAcc_Table_R;
   ExtAcc_NBin = ExtAcc_AuxStruct.ExtAcc_NBin;  

   Acc[0] = Mis_InterpolateFromTable( ExtAcc_NBin, Table_R, Table_Acc, r ) * dx / r;
   Acc[1] = Mis_InterpolateFromTable( ExtAcc_NBin, Table_R, Table_Acc, r ) * dy / r;
   Acc[2] = Mis_InterpolateFromTable( ExtAcc_NBin, Table_R, Table_Acc, r ) * dz / r;

} // FUNCTION : ExternalAcc



#endif // #ifdef GRAVITY



#endif // #ifndef __EXTERNALACC__
