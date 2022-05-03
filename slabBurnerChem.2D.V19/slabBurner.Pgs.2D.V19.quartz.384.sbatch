#!/bin/bash
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 384
#SBATCH -J 3dSlab
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb

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

export TITLE=2Dn384_G02
export VELOCITY=.212175

# export TITLE=2Dn384_G20
# export VELOCITY="t<0.06?1.0:16.967126,0.0"
# 
# export TITLE=2Dn384_G50
# export VELOCITY="t<0.06?1.0:42.417815,0.0"
# 
# export TITLE=2Dn384_G100
# export VELOCITY="t<0.06?1.0:84.835630,0.0"

##### Launch parallel job using srun
srun -n13824 /g/g15/mcgurn4/ablateOpt/ablate \
	--input /p/lustre1/mcgurn4/ablateInputs/slabBurnerChem.2D.V19/slabBurner.Pgs.2D.V19.yaml \
	-yaml::environment::title $TITLE \
	-yaml::solvers::[1]::processes::[0]::velocity \"$VELOCITY\"
	
echo 'Done'
