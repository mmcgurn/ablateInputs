#!/bin/bash
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 384
#SBATCH -J 2dSlab
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb
#SBATCH --mail-user=mtmcgurn@buffalo.edu

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

export TITLE=v0-280x40x40-384
export FACES=280,40,40
export UPWIND=upwind

# export TITLE=v0-280x40x40-384-ls
# export FACES=280,40,40
# export UPWIND=leastsquares

# export TITLE=v0-420x60x60-384
# export FACES=420,60,60
# export UPWIND=upwind

# export TITLE=v0-280x40x40-384-ls
# export FACES=420,60,60
# export UPWIND=leastsquares

##### Launch parallel job using srun
srun -n13824 /g/g15/mcgurn4/ablateOpt/ablate \
	--input /p/lustre2/mcgurn4/ablateInputs/slabBurner3D/slabBurner.3D.V0.yaml \
	-yaml::environment::title $TITLE \
	-yaml::timestepper::domain::faces [$FACES] \
	-yaml::timestepper::domain::fields::[0]::conservedFieldOptions::petscfv_type $UPWIND

	
echo 'Done'
