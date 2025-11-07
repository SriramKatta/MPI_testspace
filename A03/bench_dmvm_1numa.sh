#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=1:30:00
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

MAT_SIZES=("1000 300" "4000 300" "10000 5" "20000 2")
RES_FILE="results/intra_num_${SLURM_JOB_ID}"

echo "# iterations, procs, problem size, flop rate, walltime" > $RES_FILE

for MAT_SIZE in "${MAT_SIZES[@]}"
do
    for NP in {1..17}
    do
        echo "start for NP : $NP for mat size and iter $MAT_SIZE" 
        likwid-mpirun -mpi slurm \
                      -n $NP  ./exe-ICX $MAT_SIZE \
                      | grep "RES" \
                      | awk '{print $3, $4, $5, $6}' \
                      >> $RES_FILE
    done
done


