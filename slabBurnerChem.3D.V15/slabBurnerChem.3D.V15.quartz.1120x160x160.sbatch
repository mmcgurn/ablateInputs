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

##### Launch parallel job using srun
srun -n13824 /g/g15/mcgurn4/ablateOpt/ablate \
	--input /p/lustre1/mcgurn4/ablateInputs/slabBurnerChem.3D.V15/slabBurnerChem.3D.V15.yaml \
	-yaml::environment::title pgs_1120x160x160_upwind \
	-yaml::timestepper::domain::fields::[0]::conservedFieldOptions::petscfv_type upwind \
	-yaml::timestepper::domain::faces [1120,160,160]
	
echo 'Done'
