#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=0:05:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV

module load intel 
module load intelmpi

make clean 
make -j

srun -n 1 ./exe-ICX