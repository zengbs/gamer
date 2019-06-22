#!/bin/bash
gcc main.cpp  CPU_Shared_FluUtility.cpp \
-DCONSERVED_ENERGY=2 \
-DSERIAL \
-DFLOAT8 \
-DLR_SCHEME=PLM \
-DMODEL=SR_HYDRO \
-DRANDOM_NUMBER=RNG_GNU_EXT \
-DFLU_SCHEME=MHM \
-DRSOLVER=HLLC \
-DNLEVEL=10 \
-DMAX_PATCH=200000 \
-DEOS=APPROXIMATED_GENERAL    -lm && ./a.out
#-DEOS=CONSTANT_GAMMA         -lm && ./a.out
#-DEOS=RELATIVISTIC_IDEAL_GAS -lm && ./a.out
