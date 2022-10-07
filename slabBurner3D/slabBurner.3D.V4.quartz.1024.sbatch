#!/bin/bash
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 1024
#SBATCH -J slabV3
#SBATCH -t 64:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb
#SBATCH --mail-user=mtmcgurn@buffalo.edu
#SBATCH --reservation=dat_10_7
#SBATCH --qos=exempt

##### Load Required modules
# gcc
module load mkl/2019.0
module load valgrind/3.16.1
module load gcc/10.2.1
module load  cmake/3.21.1 

# Load PETSC ENV
export PETSC_DIR="/g/g15/mcgurn4/petsc"
export PETSC_ARCH="arch-ablate-opt-gcc" # arch-ablate-debug or arch-ablate-opt
export PKG_CONFIG_PATH="${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$PKG_CONFIG_PATH"
export HDF5_ROOT="${PETSC_DIR}/${PETSC_ARCH}"  
# Include the bin directory to access mpi commands
export PATH="${PETSC_DIR}/${PETSC_ARCH}/bin:$PATH"

# Make a temp directory so that tchem has a place to vomit its files
mkdir tmp_$SLURM_JOBID
cd tmp_$SLURM_JOBID

export TITLE=v4-490x70x70-1024
export FACES=490,70,70
export UPWIND=leastsquares

# export TITLE=v4-560x80x80-1024
# export FACES=560,80,80
# export UPWIND=leastsquares

##### Launch parallel job using srun
srun -n36864 /g/g15/mcgurn4/ablateOpt/ablate \
	--input /p/lustre2/mcgurn4/ablateInputs/slabBurner3D/slabBurner.3D.V4.yaml \
	-yaml::environment::title $TITLE \
	-yaml::timestepper::domain::faces [$FACES] \
	-yaml::timestepper::domain::fields::[0]::conservedFieldOptions::petscfv_type $UPWIND

	
echo 'Done'
