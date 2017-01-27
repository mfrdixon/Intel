#!/usr/bin/bash
export MIC_ENV_PREFIX=PHI
export PHI_KMP_AFFINITY=compact
#export STATES_C=off
export turbo=enabled
cores=( 1 2 4 8 16 32 64 68 )
tpcs=( 1 2 3 4 )
echo "\begin{table}"
echo "\begin{tabular}{|c|c|c|c|c|} \hline"
echo "Cores & 1 thread/core & 2 threads/core & 3 threads/core & 4 threads/core\\\\"
for core in "${cores[@]}"
do
 echo $core
 for tpc in "${tpcs[@]}"
 do
   echo "&"
   export KMP_HW_SUBSET=${core}c,${tpc}t
   export PHI_OMP_NUM_THREADS=$(($cores * $tpc))
   ./HestonCalibration.MIC
  done
 echo "\\\\"
done
echo "\hline"
echo "\end{tabular}"
echo "\end{table}"
