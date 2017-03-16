#!/bin/bash

# set the number of runs
RUNS=2000

# Set the default OpenMP behaviour
export OMP_WAIT_POLICY=PASSIVE

# Clear up the memory first
./memhog `free -g | grep Mem | awk '{print int(0.9*$2)}'`

        ./jac_noverlet 5dhfr_cube.psf 5dhfr_cube.pdb 16 $RUNS > jac_hybrid_noverlet_1_16.dump
        grep -v ":" jac_hybrid_noverlet_1_16.dump | grep -v "#" | awk 'NR>1{ for(i=7;i<=NF;i++){ sums[i]+=$i }} END { for(i=7;i<=NF;i++){ printf "%.3f  ", sums[i]/(NR-1)}; print "" }' >> jac_hybrid.totals

        mpirun -n 2 ./jac_mpi_noverlet 5dhfr_cube.psf 5dhfr_cube.pdb 8 $RUNS > jac_hybrid_noverlet_2_8.dump
        grep -v ":" jac_hybrid_noverlet_2_8.dump | grep -v "#" | awk 'NR>1{ for(i=7;i<=NF;i++){ sums[i]+=$i }} END { for(i=7;i<=NF;i++){ printf "%.3f  ", sums[i]/(NR-1)}; print "" }' >> jac_hybrid.totals

        mpirun -n 4 ./jac_mpi_noverlet 5dhfr_cube.psf 5dhfr_cube.pdb 4 $RUNS > jac_hybrid_noverlet_4_4.dump
        grep -v ":" jac_hybrid_noverlet_4_4.dump | grep -v "#" | awk 'NR>1{ for(i=7;i<=NF;i++){ sums[i]+=$i }} END { for(i=7;i<=NF;i++){ printf "%.3f  ", sums[i]/(NR-1)}; print "" }' >> jac_hybrid.totals

        mpirun -n 8 ./jac_mpi_noverlet 5dhfr_cube.psf 5dhfr_cube.pdb 2 $RUNS > jac_hybrid_noverlet_8_2.dump
        grep -v ":" jac_hybrid_noverlet_8_2.dump | grep -v "#" | awk 'NR>1{ for(i=7;i<=NF;i++){ sums[i]+=$i }} END { for(i=7;i<=NF;i++){ printf "%.3f  ", sums[i]/(NR-1)}; print "" }' >> jac_hybrid.totals

        mpirun -n 16 ./jac_mpi_noverlet 5dhfr_cube.psf 5dhfr_cube.pdb 1 $RUNS > jac_hybrid_noverlet_16_1.dump
        grep -v ":" jac_hybrid_noverlet_16_1.dump | grep -v "#" | awk 'NR>1{ for(i=7;i<=NF;i++){ sums[i]+=$i }} END { for(i=7;i<=NF;i++){ printf "%.3f  ", sums[i]/(NR-1)}; print "" }' >> jac_hybrid.totals

