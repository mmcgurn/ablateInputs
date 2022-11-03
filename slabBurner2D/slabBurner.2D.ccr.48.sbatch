#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --qos=supporters
#SBATCH --account=chrest
#SBATCH --time=72:00:00
#SBATCH --nodes=48
#SBATCH --constraint=IB
#SBATCH --ntasks-per-node=40
#SBATCH --job-name="slbBurn"
#SBATCH --mail-user=mtmcgurn@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --exclusive

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


export TITLE=v3-6G-560x80-48
export FACES=560,80
export UPWIND=leastsquares
export VELOCITY="min(3.985120454,t*3.985120454/.01),0.0"

# export TITLE=v3-22G-560x80-48
# export FACES=560,80
# export UPWIND=leastsquares
# export VELOCITY="min(14.96274499,t*14.96274499/.01),0.0"


# Run all tests
echo "Start Time " `date +%s`

srun -n $SLURM_NPROCS  /projects/academic/chrest/mtmcgurn/ablateOpt/ablate \
  --input /panasas/scratch/grp-chrest/mtmcgurn/ablateInputs/slabBurner2D/slabBurner.2D.yaml \
  -yaml::environment::title $TITLE \
  -yaml::timestepper::domain::faces [$FACES] \
  -yaml::timestepper::domain::fields::[0]::conservedFieldOptions::petscfv_type $UPWIND \
  -yaml::solvers::[1]::processes::[0]::velocity \"$VELOCITY\"
 

echo "End Time " `date +%s`
