#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=10:00:00
#SBATCH --cpu-freq=2000000-2000000:performance
#SBATCH --job-name=single_numa_dmvm
#SBATCH --output=./SLURM_OUT_FILES/out/%j_%x.out
#SBATCH --error=./SLURM_OUT_FILES/err/%j_%x.err
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV

module load intel
module load intelmpi
module load likwid

make

mkdir -p results

MAT_SIZES=("1000 1000000" "4000 100000" "10000 10000" "20000 5000")
RES_FILE="results/inter_numa_${SLURM_JOB_ID}.csv"

echo "# iterations, procs, problem size, flop rate, walltime" > $RES_FILE

for MAT_SIZE in "${MAT_SIZES[@]}"
do
    for NP in {18..72..18}
    do
        echo "start for NP : $NP for mat size and iter $MAT_SIZE" 
        likwid-mpirun -mpi slurm \
                      -n $NP  ./exe-ICX $MAT_SIZE \
                      | grep "RES" \
                      | awk '{printf "%s,%s,%s,%s,%s\n", $3, $4, $5, $6, $7}' \
                      >> $RES_FILE
    done
done


