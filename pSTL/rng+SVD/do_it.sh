#! /bin/bash
cores=( 1 2 4 8 16 32 68 )
threads=( 1 2 3 4)

for i in "${cores[@]}"
do
    for j in "${threads[@]}"
    do
        build/src/svd 512 $i $j
    done
done
