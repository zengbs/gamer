#include "GAMER.h"
#include "CUFLU.h"

#if ( !defined GPU  &&  MODEL == SR_HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP) )


extern real CPU_CheckMinPres( const real InPres, const real MinPres );

static void LimitSlope( const real L2[], const real L1[], const real C0[], const real R1[], const real R2[],
                        const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                        const real Gamma, const int XYZ, real Slope_Limiter[] );


#if ( LR_SCHEME == PLM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-linear method (PLM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. The input and output data should be primitive variables
//                3. The PLM and PPM data reconstruction functions share the same function name
//                4. The face-centered variables will be advanced by half time-step for the CTU scheme
//                5. The data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION"
//                6. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  PriVar         : Array storing the input primitive variables
//                FC_Var         : Array to store the output face-centered primitive variables
//                NIn            : Size of the input array "PriVar" in one direction
//                NGhost         : Size of the ghost zone
//                                  --> "NIn-2*NGhost" cells will be computed along each direction
//                                  --> The size of the output array "FC_Var" is assumed to be "(NIn-2*NGhost)^3"
//                                  --> The reconstructed data at cell (i,j,k) will be stored in the
//                                      array "FC_Var" with the index "(i-NGhost,j-NGhost,k-NGhost)
//                Gamma          : Ratio of specific heats
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//                dt             : Time interval to advance solution (for the CTU scheme)
//                dh             : Grid size (for the CTU scheme)
//                MinDens/Pres   : Minimum allowed density and pressure
//------------------------------------------------------------------------------------------------------
void CPU_DataReconstruction( const real PriVar[][NCOMP_TOTAL], real FC_Var[][6][NCOMP_TOTAL], const int NIn, const int NGhost,
                             const real Gamma, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                             const real EP_Coeff, const real dt, const real dh, const real MinDens, const real MinPres )
{

   const int dr1[3] = { 1, NIn, NIn*NIn };
   const int NOut   = NIn - 2*NGhost;                    // number of output grids
   int  ID1, ID2, ID1_L, ID1_R, ID1_LL, ID1_RR, dL, dR;
   real Min, Max;
   real Slope_Limiter[NCOMP_TOTAL] = { (real)0.0 };


   for (int k1=NGhost, k2=0;  k1<NGhost+NOut;  k1++, k2++)
   for (int j1=NGhost, j2=0;  j1<NGhost+NOut;  j1++, j2++)
   for (int i1=NGhost, i2=0;  i1<NGhost+NOut;  i1++, i2++)
   {
      ID1 = (k1*NIn  + j1)*NIn  + i1;
      ID2 = (k2*NOut + j2)*NOut + i2;

//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {

//       (2-1) evaluate the monotonic slope
         dL    = 2*d;
         dR    = dL+1;
         ID1_L = ID1 - dr1[d];
         ID1_R = ID1 + dr1[d];

         if ( LR_Limiter == EXTPRE )
         {
            ID1_LL = ID1 - 2*dr1[d];
            ID1_RR = ID1 + 2*dr1[d];

            LimitSlope( PriVar[ID1_LL], PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], PriVar[ID1_RR], LR_Limiter,
                        MinMod_Coeff, EP_Coeff, Gamma, d, Slope_Limiter );
         }

         else
         {
            LimitSlope( NULL, PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], NULL, LR_Limiter,
                        MinMod_Coeff, NULL_REAL, Gamma, d, Slope_Limiter );
         }


//       (2-2) get the face-centered primitive variables
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            FC_Var[ID2][dL][v] = PriVar[ID1][v] - (real)0.5*Slope_Limiter[v];
            FC_Var[ID2][dR][v] = PriVar[ID1][v] + (real)0.5*Slope_Limiter[v];
         }


//       (2-3) ensure the face-centered variables lie between neighboring cell-centered values
         if ( LR_Limiter != EXTPRE )
         {
            for (int v=0; v<NCOMP_TOTAL; v++)
            {
               Min = ( PriVar[ID1][v] < PriVar[ID1_L][v] ) ? PriVar[ID1][v] : PriVar[ID1_L][v];
               Max = ( PriVar[ID1][v] > PriVar[ID1_L][v] ) ? PriVar[ID1][v] : PriVar[ID1_L][v];
               FC_Var[ID2][dL][v] = ( FC_Var[ID2][dL][v] > Min  ) ? FC_Var[ID2][dL][v] : Min;
               FC_Var[ID2][dL][v] = ( FC_Var[ID2][dL][v] < Max  ) ? FC_Var[ID2][dL][v] : Max;
               FC_Var[ID2][dR][v] = (real)2.0*PriVar[ID1][v] - FC_Var[ID2][dL][v];

               Min = ( PriVar[ID1][v] < PriVar[ID1_R][v] ) ? PriVar[ID1][v] : PriVar[ID1_R][v];
               Max = ( PriVar[ID1][v] > PriVar[ID1_R][v] ) ? PriVar[ID1][v] : PriVar[ID1_R][v];
               FC_Var[ID2][dR][v] = ( FC_Var[ID2][dR][v] > Min  ) ? FC_Var[ID2][dR][v] : Min;
               FC_Var[ID2][dR][v] = ( FC_Var[ID2][dR][v] < Max  ) ? FC_Var[ID2][dR][v] : Max;
               FC_Var[ID2][dL][v] = (real)2.0*PriVar[ID1][v] - FC_Var[ID2][dR][v];
            }
         }

         else // for the extrema-preserving limiter --> ensure positive density and pressure
         {
            FC_Var[ID2][dL][0] = FMAX( FC_Var[ID2][dL][0], MinDens );
            FC_Var[ID2][dR][0] = FMAX( FC_Var[ID2][dR][0], MinDens );

            FC_Var[ID2][dL][4] = CPU_CheckMinPres( FC_Var[ID2][dL][4], MinPres );
            FC_Var[ID2][dR][4] = CPU_CheckMinPres( FC_Var[ID2][dR][4], MinPres );
         }

      } // for (int d=0; d<3; d++)
   } // k,j,i

} // FUNCTION : CPU_DataReconstruction (PLM)
#endif // #if ( LR_SCHEME == PLM )



#if ( LR_SCHEME == PPM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-parabolic method (PPM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. The input and output data should be primitive variables
//                3. The PLM and PPM data reconstruction functions share the same function name
//                4. The face-centered variables will be advanced by half time-step for the CTU scheme
//                5. Currently the extrema-preserving limiter is not supported in PPM
//                6. The data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION"
//                7. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  PriVar         : Array storing the input primitive variables
//                FC_Var         : Array to store the output face-centered primitive variables
//                NIn            : Size of the input array "PriVar" in one direction
//                NGhost         : Size of the ghost zone
//                                  --> "NIn-2*NGhost" cells will be computed along each direction
//                                  --> The size of the output array "FC_Var" is assumed to be "(NIn-2*NGhost)^3"
//                                  --> The reconstructed data at cell (i,j,k) will be stored in the
//                                      array "FC_Var" with the index "(i-NGhost,j-NGhost,k-NGhost)
//                Gamma          : Ratio of specific heats
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter (useless in PPM)
//                dt             : Time interval to advance solution (for the CTU scheme)
//                dh             : Grid size (for the CTU scheme)
//                MinDens/Pres   : Minimum allowed density and pressure
//------------------------------------------------------------------------------------------------------
void CPU_DataReconstruction( const real PriVar[][NCOMP_TOTAL], real FC_Var[][6][NCOMP_TOTAL], const int NIn, const int NGhost,
                             const real Gamma, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                             const real EP_Coeff, const real dt, const real dh, const real MinDens, const real MinPres )
{

// check
#  ifdef GAMER_DEBUG
   if ( LR_Limiter == EXTPRE )
      Aux_Error( ERROR_INFO, "PPM reconstruction does NOT support the extrema-preserving limiter !!\n");
#  endif


   const int NOut   = NIn - 2*NGhost;                    // number of output grids
   const int NSlope = NOut + 2;                          // number of grids required to store the slope data
   const int dr1[3] = { 1, NIn, NIn*NIn };
   const int dr3[3] = { 1, NSlope, NSlope*NSlope };

   int ID1, ID2, ID3, ID1_L, ID1_R, ID3_L, ID3_R, dL, dR;
   real Slope_Limiter[NCOMP_TOTAL] = { (real)0.0 };
   real CC_L, CC_R, CC_C, dCC_L, dCC_R, dCC_C, FC_L, FC_R, dFC[NCOMP_TOTAL], dFC6[NCOMP_TOTAL], Max, Min;

   real (*Slope_PPM)[3][NCOMP_TOTAL] = new real [ NSlope*NSlope*NSlope ][3][NCOMP_TOTAL];



// (2-1) evaluate the monotonic slope
   for (int k1=NGhost-1, k2=0;  k1<NGhost-1+NSlope;  k1++, k2++)
   for (int j1=NGhost-1, j2=0;  j1<NGhost-1+NSlope;  j1++, j2++)
   for (int i1=NGhost-1, i2=0;  i1<NGhost-1+NSlope;  i1++, i2++)
   {
      ID1 = (k1*NIn    + j1)*NIn    + i1;
      ID2 = (k2*NSlope + j2)*NSlope + i2;

//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         ID1_L = ID1 - dr1[d];
         ID1_R = ID1 + dr1[d];

         if ( LR_Limiter == EXTPRE )
         {
            Aux_Error( ERROR_INFO, "PPM reconstruction does NOT support the extrema-preserving limiter !!\n");

            /*
            ID1_LL = ID1 - 2*dr1[d];
            ID1_RR = ID1 + 2*dr1[d];

            LimitSlope( PriVar[ID1_LL], PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], PriVar[ID1_RR], LR_Limiter,
                        MinMod_Coeff, EP_Coeff, Gamma, d, Slope_Limiter );
            */
         }

         else
         {
            LimitSlope( NULL, PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], NULL, LR_Limiter,
                        MinMod_Coeff, NULL_REAL, Gamma, d, Slope_Limiter );
         }


//       store the slope to the array "Slope_PPM"
         for (int v=0; v<NCOMP_TOTAL; v++)   Slope_PPM[ID2][d][v] = Slope_Limiter[v];

      } // for (int d=0; d<3; d++)
   } // k,j,i


   for (int k1=NGhost, k2=0, k3=1;  k1<NGhost+NOut;  k1++, k2++, k3++)
   for (int j1=NGhost, j2=0, j3=1;  j1<NGhost+NOut;  j1++, j2++, j3++)
   for (int i1=NGhost, i2=0, i3=1;  i1<NGhost+NOut;  i1++, i2++, i3++)
   {
      ID1 = (k1*NIn    + j1)*NIn    + i1;
      ID2 = (k2*NOut   + j2)*NOut   + i2;
      ID3 = (k3*NSlope + j3)*NSlope + i3;


//    (2-3) get the face-centered primitive variables
//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         dL    = 2*d;
         dR    = dL+1;
         ID1_L = ID1 - dr1[d];
         ID1_R = ID1 + dr1[d];
         ID3_L = ID3 - dr3[d];
         ID3_R = ID3 + dr3[d];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
//          (2-3-1) parabolic interpolation
            CC_L  = PriVar[ID1_L][v];
            CC_R  = PriVar[ID1_R][v];
            CC_C  = PriVar[ID1  ][v];

            dCC_L = Slope_PPM[ID3_L][d][v];
            dCC_R = Slope_PPM[ID3_R][d][v];
            dCC_C = Slope_PPM[ID3  ][d][v];

            FC_L  = (real)0.5*( CC_C + CC_L ) - (real)1.0/(real)6.0*( dCC_C - dCC_L );
            FC_R  = (real)0.5*( CC_C + CC_R ) - (real)1.0/(real)6.0*( dCC_R - dCC_C );


//          (2-3-2) monotonicity constraint
            dFC [v] = FC_R - FC_L;
            dFC6[v] = (real)6.0*(  CC_C - (real)0.5*( FC_L + FC_R )  );

            if (  ( FC_R - CC_C )*( CC_C - FC_L ) <= (real)0.0  )
            {
               FC_L = CC_C;
               FC_R = CC_C;
            }
            else if ( dFC[v]*dFC6[v] > +dFC[v]*dFC[v] )
               FC_L = (real)3.0*CC_C - (real)2.0*FC_R;
            else if ( dFC[v]*dFC6[v] < -dFC[v]*dFC[v] )
               FC_R = (real)3.0*CC_C - (real)2.0*FC_L;


//          (2-3-3) ensure the face-centered variables lie between neighboring cell-centered values
            Min  = ( CC_C < CC_L ) ? CC_C : CC_L;
            Max  = ( CC_C > CC_L ) ? CC_C : CC_L;
            FC_L = ( FC_L > Min  ) ? FC_L : Min;
            FC_L = ( FC_L < Max  ) ? FC_L : Max;

            Min  = ( CC_C < CC_R ) ? CC_C : CC_R;
            Max  = ( CC_C > CC_R ) ? CC_C : CC_R;
            FC_R = ( FC_R > Min  ) ? FC_R : Min;
            FC_R = ( FC_R < Max  ) ? FC_R : Max;


            FC_Var[ID2][dL][v] = FC_L;
            FC_Var[ID2][dR][v] = FC_R;

         } // for (int v=0; v<NCOMP_TOTAL; v++)

      } // for (int d=0; d<3; d++)
   } // k,j,i

   delete [] Slope_PPM;

} // FUNCTION : CPU_DataReconstruction (PPM)
#endif // #if ( LR_SCHEME == PPM )




//-------------------------------------------------------------------------------------------------------
// Function    :  LimitSlope
// Description :  Evaluate the monotonic slope by applying slope limiters
//
// Note        :  1. The input data should be primitive variables
//                2. The L2 and R2 elements are useful only for the extrema-preserving limiter
//
// Parameter   :  L2             : Element x-2
//                L1             : Element x-1
//                C0             : Element x
//                R1             : Element x+1
//                R2             : Element x+2
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//                Gamma          : Ratio of specific heats
//                                 --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                XYZ            : Target spatial direction : (0/1/2) --> (x/y/z)
//                                 --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                Slope_Limiter  : Array to store the output monotonic slope
//-------------------------------------------------------------------------------------------------------
void LimitSlope( const real L2[], const real L1[], const real C0[], const real R1[], const real R2[],
                 const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                 const real Gamma, const int XYZ, real Slope_Limiter[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( LR_Limiter == EXTPRE  &&  ( L2 == NULL || R2 == NULL )  )
      Aux_Error( ERROR_INFO, "input element == NULL !!\n" );
#  endif


   real Slope_L[NCOMP_TOTAL], Slope_R[NCOMP_TOTAL], Slope_C[NCOMP_TOTAL], Slope_A[NCOMP_TOTAL];
   real Slope_LL[NCOMP_TOTAL], Slope_RR[NCOMP_TOTAL], Slope_LR;
   real D2_L, D2_R, D2_C, D2_Sign, D2_Limiter, Slope_Sign;  // variables for the extrema-preserving limiter


// evaluate different slopes
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      Slope_L[v] = C0[v] - L1[v];
      Slope_R[v] = R1[v] - C0[v];
      Slope_C[v] = (real)0.5*( Slope_L[v] + Slope_R[v] );
   }

   if ( LR_Limiter == VL_GMINMOD )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         if ( Slope_L[v]*Slope_R[v] > (real)0.0 )
            Slope_A[v] = (real)2.0*Slope_L[v]*Slope_R[v]/( Slope_L[v] + Slope_R[v] );
         else
            Slope_A[v] = (real)0.0;
      }
   }

   if ( LR_Limiter == EXTPRE )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         Slope_LL[v] = L1[v] - L2[v];
         Slope_RR[v] = R2[v] - R1[v];
      }
   }


// apply the slope limiter
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      Slope_LR = Slope_L[v]*Slope_R[v];

      if (  Slope_LR > (real)0.0  &&  ( LR_Limiter != EXTPRE || Slope_LL[v]*Slope_RR[v] > (real)0.0 )  )
      {
         switch ( LR_Limiter )
         {
            case VANLEER:              // van-Leer
               Slope_Limiter[v] = (real)2.0*Slope_LR/( Slope_L[v] + Slope_R[v] );
               break;

            case GMINMOD: case EXTPRE: // generalized MinMod & extrema-preserving
               Slope_L[v] *= MinMod_Coeff;
               Slope_R[v] *= MinMod_Coeff;
               Slope_Limiter[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_C[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v] *= SIGN( Slope_C[v] );
               break;

            case ALBADA:               // van-Albada
               Slope_Limiter[v] = Slope_LR*( Slope_L[v] + Slope_R[v] ) /
                                  ( Slope_L[v]*Slope_L[v] + Slope_R[v]*Slope_R[v] );
               break;

            case VL_GMINMOD:           // van-Leer + generalized MinMod
               Slope_L[v] *= MinMod_Coeff;
               Slope_R[v] *= MinMod_Coeff;
               Slope_Limiter[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_C[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_A[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v] *= SIGN( Slope_C[v] );
               break;

            default :
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "LR_Limiter", LR_Limiter );
         }
      } // if (  Slope_LR > (real)0.0  &&  ( LR_Limiter != EXTPRE || Slope_LL[v]*Slope_RR[v] > (real)0.0 )  )

      else
      {
         if ( LR_Limiter == EXTPRE )   // extrema-preserving
         {
            D2_L = Slope_L [v] - Slope_LL[v];
            D2_R = Slope_RR[v] - Slope_R [v];
            D2_C = Slope_R [v] - Slope_L [v];

            D2_Sign    = SIGN( D2_C );
            Slope_Sign = SIGN( Slope_C[v] );

            D2_Limiter = FMIN(  FABS(D2_C), FMIN( FMAX(D2_Sign*D2_L, (real)0.0),
                                                  FMAX(D2_Sign*D2_R, (real)0.0) )  );

            if ( D2_Sign*Slope_Sign < (real)0.0 )
               Slope_Limiter[v] = FMIN( (real)1.5*EP_Coeff*D2_Limiter, MinMod_Coeff*FABS(Slope_L[v]) );
            else
               Slope_Limiter[v] = FMIN( (real)1.5*EP_Coeff*D2_Limiter, MinMod_Coeff*FABS(Slope_R[v]) );

            Slope_Limiter[v] = Slope_Sign * FMIN( FABS(Slope_C[v]), Slope_Limiter[v] );
         }
         else
            Slope_Limiter[v] = (real)0.0;

      } // if ( Slope_LR > (real)0.0 && ( LR_Limiter != EXTPRE || Slope_LL[v]*Slope_RR[v] > (real)0.0 ) ) .else.
   } // for (int v=0; v<NCOMP_TOTAL; v++)

} // FUNCTION : LimitSlope



#endif // #if ( !defined GPU  &&  MODEL == SR_HYDRO  &&  (FLU_SCHEME == MHM || MHM_RP) )