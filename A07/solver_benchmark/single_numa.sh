#!/bin/bash -l
#
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH -J single_numa_scale
#SBATCH --output=SLURM_OUT_FILES/%x_%j.out
#SBATCH --error=SLURM_ERR_FILES/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV


module load intel
module load intelmpi
module load likwid

DIR_NAME=benchmark_${SLURM_JOB_NAME}_${SLURM_JOB_ID}
EXE_NAME=EXE_${SLURM_JOB_ID}

mkdir -p $DIR_NAME
cp ./exe-ICX ./${DIR_NAME}/${EXE_NAME}

cd $DIR_NAME

for numproc in {1..18}
do
    mkdir -p ${numproc}_cores
    echo "****************************************************************************************"
    echo "running on ${numproc} core of ${SLURM_JOB_NAME}"
    echo "****************************************************************************************"
    cd ${numproc}_cores
    likwid-mpirun -mpi slurm \
                    -omp intel \
                    -n ${numproc} \
                    ../${EXE_NAME} ../../canal_memdomain.par | tee resfile
    cd ..
done

rm ${EXE_NAME}

