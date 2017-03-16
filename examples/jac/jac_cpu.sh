#!/bin/bash

# set the number of runs
RUNS=2000

# Set the default OpenMP behaviour
export OMP_WAIT_POLICY=PASSIVE

# Set the path for CUDA libraries
export LD_LIBRARY_PATH=/opt/cuda/lib64:/usr/local/cuda/lib64:/usr/local/cuda/lib:$LD_LIBRARY_PATH

# loop over number of CPUs
for cpu in {1..12}
do

    # Set some environment variables for OpenMP
    export OMP_NUM_THREADS=$cpu
    export OMP_THREAD_LIMIT=$cpu
    export OMP_PROC_BIND=TRUE
    
    # The threaded runs on a single node
    for vtype in noverlet verlet pwverlet pwverlet2
    do
    
        ./jac_${vtype} 5dhfr_cube.psf 5dhfr_cube.pdb $cpu $RUNS > jac_${vtype}_${cpu}.dump
        grep -v ":" jac_${vtype}_${cpu}.dump | grep -v "#" | awk 'NR>1{ for(i=7;i<=NF;i++){ sums[i]+=$i }} END { for(i=7;i<=NF;i++){ printf "%.3f  ", sums[i]/(NR-1)}; print "" }' >> jac_${vtype}.totals

    done

done
