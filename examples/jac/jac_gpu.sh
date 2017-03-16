#!/bin/bash

# set the number of runs
RUNS=2000

# Set the default OpenMP behaviour
export OMP_WAIT_POLICY=PASSIVE
export OMP_NUM_THREADS=2

# Set the path for CUDA libraries
export LD_LIBRARY_PATH=/opt/cuda/lib64:/usr/local/cuda/lib64:/usr/local/cuda/lib:$LD_LIBRARY_PATH

# loop over number of cores
for cores in {1..128}
do

    # loop over the verlet types
    for vtype in noverlet verlet
    do

        # The GPU runs on a single node
        ./jac_cuda_${vtype} 5dhfr_cube.psf 5dhfr_cube.pdb $cores $RUNS > jac_cuda_${vtype}_${cores}.dump
        grep -v ":" jac_cuda_${vtype}_${cores}.dump | grep -v "#" | awk 'NR>1{ for(i=7;i<=NF;i++) { sums[i]+=$i }} END { for(i=7;i<=NF;i++){ printf "%.3f  ", sums[i]/(NR-1)}; print "" }' >> jac_cuda_${vtype}.totals
        grep timers jac_cuda_${vtype}_${cores}.dump | awk 'NR>1{ for(i=0;i<=NF;i++) { sums[i]+=$i }} END { for(i=0;i<=NF;i++){ printf "%.3f  ", sums[i]/(NR-1)}; print "" }' >> jac_cuda_${vtype}.timers

    done

done
