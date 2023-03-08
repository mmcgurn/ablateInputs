#!/bin/bash
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 24
#SBATCH -J 6G-2D
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

export TITLE=6G-186x40-pmma-rad-soot-$SLURM_JOBID
export FACES=186,40
export FILE=/p/lustre2/ubchrest/ablateInputs/slabBurner.pmma.2D/slabBurner.2D.6G.pmma.rad.soot.yaml

# export TITLE=6G-280x60-pmma-rad-soot-$SLURM_JOBID
# export FACES=280,60
# export FILE=/p/lustre2/ubchrest/ablateInputs/slabBurner.pmma.2D/slabBurner.2D.6G.pmma.rad.soot.yaml
# 
# export TITLE=6G-372x80-pmma-rad-soot-$SLURM_JOBID
# export FACES=372,80
# export FILE=/p/lustre2/ubchrest/ablateInputs/slabBurner.pmma.2D/slabBurner.2D.6G.pmma.rad.soot.yaml


##### Launch parallel job using srun
srun -n1344 /p/lustre2/mcgurn4/ablateOpt/ablate \
   --input $FILE \
   -yaml::environment::title $TITLE \
   -yaml::timestepper::domain::faces [$FACES]

echo 'Done'
