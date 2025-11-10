#!/bin/bash -l
#
#SBATCH --nodes=4
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

export I_MPI_PIN=1
export I_MPI_DEBUG=0

NPM=18

make

mkdir -p results

MAT_SIZES=("1000 1000000" "4000 100000" "10000 10000" "20000 5000")
RES_FILE="results/inter_numa_${SLURM_JOB_ID}.csv"

echo "# iterations, procs, problem size, flop rate, walltime" > $RES_FILE

for MAT_SIZE in "${MAT_SIZES[@]}"
do
    for NP in {1..4}
    do
        npn=$(($np * $NPM))
        np_1=$(($npn - 1))
        export I_MPI_PIN_PROCESSOR_LIST=0-$np_1
        echo "start for NP : $NP for mat size and iter $MAT_SIZE" 
        mpirun -n $npn  ./exe-ICX $MAT_SIZE \
              | grep "RES" \
              | awk '{printf "%s,%s,%s,%s,%s\n", $3, $4, $5, $6, $7}' \
              >> $RES_FILE
    done
done


