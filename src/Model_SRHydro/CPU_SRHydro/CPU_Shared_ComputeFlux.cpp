#include "GAMER.h"
#include "CUFLU.h"
#include "../../../include/CPU_prototypes.h"

#if ( !defined GPU  &&  MODEL == SR_HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP) )

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_ComputeFlux
// Description :  Compute the face-centered fluxes by Riemann solver
//
// Note        :  1. Currently support the exact, HLLC, HLLE, and Roe solvers
//                2. The size of the input array "FC_Var" is assumed to be N_FC_VAR^3
//                   --> "N_FC_VAR-1" fluxes will be computed along each direction
//
// Parameter   : [ 1] FC_Var          : Array storing the input face-centered conserved variables
//               [ 2] FC_Flux         : Array to store the output face-centered flux
//               [ 3] NFlux           : Size of the array FC_Flux in each direction (must be >= N_FC_VAR-1)
//                                      --> The (i,j,k) flux will be stored in the array FC_Flux with
//                                          the index "(k*NFlux+j)*NFlux+i"
//                                      --> The (i,j,k) FC_Flux_x/y/z are defined at the +x/+y/+z surfaces of the
//                                          cell (i,j,k)
//               [ 4] Gap             : Number of grids to be skipped in the transverse direction
//                                      --> "(N_FC_VAR-2*Gap)^2" fluxes will be computed on each surface
//               [ 5] Gamma           : Ratio of specific heats
//               [ 6] CorrHalfVel     : true --> correcting the half-step velocity by gravity (for UNSPLIT_GRAVITY only)
//               [ 7] Pot_USG         : Array storing the input potential for CorrHalfVel     (for UNSPLIT_GRAVITY only)
//                                      --> must have the same size as FC_Var ( (PS2+2)^3 )
//               [ 8] Corner          : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//               [ 9] dt              : Time interval to advance the full-step solution       (for UNSPLIT_GRAVITY only)
//               [10] dh              : Grid size                                             (for UNSPLIT_GRAVITY only)
//               [11] Time            : Current physical time                                 (for UNSPLIT_GRAVITY only)
//               [12] GravityType     : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//               [13] ExtAcc_AuxArray : Auxiliary array for adding external acceleration          (for UNSPLIT_GRAVITY only)
//               [14] MinPres         : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void
CPU_ComputeFlux(const real FC_Var[][6][NCOMP_TOTAL],
		real FC_Flux[][3][NCOMP_TOTAL],
		const int NFlux,
		const int Gap,
		const real Gamma,
		const bool CorrHalfVel,
		const real Pot_USG[],
		const double Corner[],
		const real dt, const real dh, const double Time,
		const OptGravityType_t GravityType, const double ExtAcc_AuxArray[], const real MinPres)
{
// check
#ifdef GAMER_DEBUG
    if (CorrHalfVel)
	Aux_Error(ERROR_INFO, "CorrHalfVel is NOT supported when UNSPLIT_GRAVITY is off !!\n");
#endif				// #ifdef GAMER_DEBUG

    const int dID2[3] = { 1, N_FC_VAR, N_FC_VAR * N_FC_VAR };

    real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL];

/*----- ID1: index of flux
 *----- ID2: index of face-centered variables ------*/

    int ID1, ID2, dL, dR, start2[3] = { 0 }, end1[3] = { 0 };

// loop over different spatial directions
    for (int d = 0; d < 3; d++) {
	dL = 2 * d;
	dR = dL + 1;

	switch (d) {
	case 0:
	    start2[0] = 0;
	    start2[1] = Gap;
	    start2[2] = Gap;
	    end1[0] = N_FC_VAR - 1;
	    end1[1] = N_FC_VAR - 2 * Gap;
	    end1[2] = N_FC_VAR - 2 * Gap;
	    break;

	case 1:
	    start2[0] = Gap;
	    start2[1] = 0;
	    start2[2] = Gap;
	    end1[0] = N_FC_VAR - 2 * Gap;
	    end1[1] = N_FC_VAR - 1;
	    end1[2] = N_FC_VAR - 2 * Gap;
	    break;

	case 2:
	    start2[0] = Gap;
	    start2[1] = Gap;
	    start2[2] = 0;
	    end1[0] = N_FC_VAR - 2 * Gap;
	    end1[1] = N_FC_VAR - 2 * Gap;
	    end1[2] = N_FC_VAR - 1;
	    break;
	}

	for (int k1 = 0, k2 = start2[2]; k1 < end1[2]; k1++, k2++)
	    for (int j1 = 0, j2 = start2[1]; j1 < end1[1]; j1++, j2++)
		for (int i1 = 0, i2 = start2[0]; i1 < end1[0]; i1++, i2++) {
		    ID1 = (k1 * NFlux + j1) * NFlux + i1; // index of 
		    ID2 = (k2 * N_FC_VAR + j2) * N_FC_VAR + i2;

		    for (int v = 0; v < NCOMP_TOTAL; v++) {
			ConVar_L[v] = FC_Var[ID2][dR][v];
			ConVar_R[v] = FC_Var[ID2 + dID2[d]][dL][v];
		    }

#                   ifdef CHECK_NEGATIVE_IN_FLUID
		    if(CPU_CheckUnphysical(ConVar_L, NULL)) Aux_Message(stderr,"\nUnphysical varibles!\nfunction: %s: %d\n", __FUNCTION__, __LINE__);
		    if(CPU_CheckUnphysical(ConVar_R, NULL)) Aux_Message(stderr,"\nUnphysical varibles!\nfunction: %s: %d\n", __FUNCTION__, __LINE__);
#                   endif

#                   if ( RSOLVER == HLLC )
		    CPU_RiemannSolver_HLLC(d, FC_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres);
#                   elif ( RSOLVER == HLLE )
		    CPU_RiemannSolver_HLLE(d, FC_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres);
#                   else
#                   error : ERROR : unsupported Riemann solver !!
#                   endif
		}// i,j,k
    }// for (int d=0; d<3; d++)

}// FUNCTION : CPU_ComputeFlux

#                   endif// #if ( !defined GPU  &&  MODEL == SR_HYDRO  &&  (FLU_SCHEME == MHM || MHM_RP || CTU) )