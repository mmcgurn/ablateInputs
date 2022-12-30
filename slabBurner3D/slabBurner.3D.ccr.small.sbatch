#!/bin/sh
#SBATCH --partition=scavenger
#SBATCH --qos=scavenger
#SBATCH --account=chrest
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --constraint=IB
#SBATCH --job-name="slbBurn"
#SBATCH --mail-user=mtmcgurn@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --requeue

# Print the current environment
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

# Load the required modules
module load intel/20.2
module load intel-mpi/2020.2
module load gcc/11.2.0
module load cmake/3.22.3
module load valgrind/3.14.0
module load gdb/7.8

# setup petsc
export PETSC_DIR=/projects/academic/chrest/mtmcgurn/petsc  
export PETSC_ARCH=arch-ablate-opt
# export PETSC_ARCH=arch-ablate-debug
export PKG_CONFIG_PATH="${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$PKG_CONFIG_PATH"

# The initial srun will trigger the SLURM prologue on the compute nodes.
NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo NPROCS=$NPROCS

# The PMI library is necessary for srun
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Tell the tests what mpi command to use
export TEST_MPI_COMMAND=srun

# Make a temp directory so that tchem has a place to vomit its files
mkdir tmp_$SLURM_JOBID
cd tmp_$SLURM_JOBID

export TITLE=6G-600-$SLURM_JOBID
export VELOCITY="min(3.985120454,t*3.985120454/.001),0.0,0.0"

##### Launch parallel job using srun
srun -n 56 \
 	/projects/academic/chrest/mtmcgurn/ablateOpt/ablate \
	--input /panasas/scratch/grp-chrest/mtmcgurn/ablateInputs/slabBurner3D/slabBurner.3D.yaml \
	-yaml::environment::title $TITLE \
	-yaml::timestepper::domain::faces [112,17,17] \
    -yaml::solvers::[1]::processes::[0]::velocity \"$VELOCITY\" \
	-yaml::timestepper::io::interval 100
	
#     
# srun -n 56 --overcommit \
#   /projects/academic/chrest/valgrind/bin/valgrind --tool=memcheck -q --num-callers=20 \
# 	--suppressions=$PETSC_DIR/share/petsc/valgrind/petsc-val.supp \
# 	--log-file=valgrind.log.%p \
#  	/projects/academic/chrest/mtmcgurn/ablateDebug/ablate \
# 	--input /panasas/scratch/grp-chrest/mtmcgurn/ablateInputs/slabBurner3D/slabBurner.3D.yaml \
# 	-yaml::timestepper::domain::faces [112,17,17] \
#     -yaml::timestepper::arguments::ts_max_steps 150 

