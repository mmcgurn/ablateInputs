#!/bin/bash
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 24
#SBATCH -J 6G-gMsh-2D
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
module load cmake/3.21.1 

# Load PETSC ENV
export PETSC_DIR="/p/lustre2/mcgurn4/petsc"
export PETSC_ARCH="arch-ablate-opt-gcc" # arch-ablate-debug or arch-ablate-opt
export PKG_CONFIG_PATH="${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$PKG_CONFIG_PATH"
export HDF5_ROOT="${PETSC_DIR}/${PETSC_ARCH}"  
# Include the bin directory to access mpi commands
export PATH="${PETSC_DIR}/${PETSC_ARCH}/bin:$PATH"

# Make a temp directory so that tchem has a place to vomit its files
mkdir tmp_$SLURM_JOBID
cd tmp_$SLURM_JOBID

export DM_REFINE=0
export TITLE=lowG-gMsh-dm$DM_REFINE-pmma-$SLURM_JOBID
export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.2D/slabBurner2D.lowG.3_8_23.yaml

# export DM_REFINE=0
# export TITLE=highG-gMsh-dm$DM_REFINE-pmma-$SLURM_JOBID
# export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.2D/slabBurner2D.highG.3_8_23.yaml


##### Launch parallel job using srun
srun -n864
 /p/lustre2/mcgurn4/ablateOpt/ablate \
   --input $FILE \
   -yaml::environment::title $TITLE \
   -yaml::timestepper::domain::options::dm_refine $DM_REFINE

echo 'Done'
