#include "CUAPI.h"

#if ( defined GPU  &&  defined GRAVITY )


#include "CUPOT.h"

extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern ExtAcc_AuxStruct_t ExtAcc_AuxStruct; 


// declare all GPU kernels requiring ExtPot_AuxArray[] and/or ExtAcc_AuxStruct
#if ( MODEL == HYDRO || MODEL == SR_HYDRO )
int CUPOT_SetConstMem_HydroGravitySolver( ExtAcc_AuxStruct_t ExtAcc_AuxStruct );
int CUPOT_SetConstMem_dtSolver_HydroGravity( ExtAcc_AuxStruct_t ExtAcc_AuxStruct );
#if (  defined UNSPLIT_GRAVITY  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )
int CUFLU_SetConstMem_FluidSolver_ExtAcc( ExtAcc_AuxStruct_t ExtAcc_AuxStruct );
#endif
#endif // #if ( MODEL == HYDRO )

#if ( MODEL == ELBDM )
int CUPOT_SetConstMem_ELBDMGravitySolver( double ExtPot_AuxArray_h[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Init_ExternalAccPot
// Description :  Set the auxiliary GPU constant-memory arrays for the external acceleration and potential

// Note        :  1. Invoked by Init_ExternalAccPot()
//
// Parameter   :  None
//
// Return      :  ExtPot_AuxArray[] and ExtAcc_AuxStruct in various CUDA kernels
//-------------------------------------------------------------------------------------------------------
void CUAPI_Init_ExternalAccPot()
{

#  if ( MODEL == HYDRO || MODEL == SR_HYDRO)
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   {
      if (  CUPOT_SetConstMem_HydroGravitySolver( ExtAcc_AuxStruct ) != 0  )
         Aux_Error( ERROR_INFO, "CUPOT_SetConstMem_HydroGravitySolver failed ...\n" );

      if (  CUPOT_SetConstMem_dtSolver_HydroGravity( ExtAcc_AuxStruct ) != 0  )
         Aux_Error( ERROR_INFO, "CUPOT_SetConstMem_dtSolver_HydroGravity failed ...\n" );

#     if (  defined UNSPLIT_GRAVITY  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )
      if (  CUFLU_SetConstMem_FluidSolver_ExtAcc( ExtAcc_AuxStruct ) != 0  )
         Aux_Error( ERROR_INFO, "CUFLU_SetConstMem_FluidSolver_ExtAcc failed ...\n" );
#     endif
   }
#  endif // if ( MODEL == HYDRO )

#  if ( MODEL == ELBDM )
   if (  OPT__EXTERNAL_POT  &&  CUPOT_SetConstMem_ELBDMGravitySolver( ExtPot_AuxArray ) != 0  )
      Aux_Error( ERROR_INFO, "CUPOT_SetConstMem_ELBDMGravitySolver failed ...\n" );
#  endif

} // FUNCTION : CUAPI_Init_ExternalAccPot



#endif // #if ( defined GPU  &&  defined GRAVITY )
