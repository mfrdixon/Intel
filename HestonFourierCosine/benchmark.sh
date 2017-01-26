#!/usr/bin/bash
export MIC_ENV_PREFIX=PHI
export PHI_KMP_AFFINITY=compact
export STATES_C=off
export turbo=enabled
array=( 1 2 4 8 16 32 60 )
tpc=2
for cores in "${array[@]}"
do
 export KMP_HW_SUBSET=${cores}c,${tpc}t
 echo $KMP_HW_SUBSET
 export PHI_OMP_NUM_THREADS=$(($cores * $tpc))
 echo $PHI_OMP_NUM_THREADS
 echo "Cores: $cores, threads per core: $tpc"
 ./HestonCalibration.MIC
done
