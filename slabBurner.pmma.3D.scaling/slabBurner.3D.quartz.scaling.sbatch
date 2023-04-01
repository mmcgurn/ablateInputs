#!/bin/bash
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 2
#SBATCH -J scale
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
export PETSC_DIR="/usr/workspace/mcgurn4/petsc"
export PETSC_ARCH="arch-ablate-opt-gcc" # arch-ablate-debug or arch-ablate-opt
export PKG_CONFIG_PATH="${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$PKG_CONFIG_PATH"
export HDF5_ROOT="${PETSC_DIR}/${PETSC_ARCH}"  
# Include the bin directory to access mpi commands
export PATH="${PETSC_DIR}/${PETSC_ARCH}/bin:$PATH"

# Make a temp directory so that tchem has a place to vomit its files
mkdir tmp_$SLURM_JOBID
cd tmp_$SLURM_JOBID

export FILE=/p/lustre2/ubchrest/ablateInputs/slabBurner.pmma.3D.scaling/slabBurner.3D.6G.pmma.rad.yaml

# 2 nodes
spindle --location=/var/tmp/mcgurn4 srun -n72 /usr/workspace/mcgurn4/ablateOpt/ablate \
   --input $FILE \
   -yaml::environment::title dm0-$SLURM_JOBID \
   -yaml::timestepper::domain::faces [140,30,30] \
   -log_view :/p/lustre2/ubchrest/ablateInputsslabBurner.pmma.3D.scaling/ConvEst_Refinement_${SLURM_JOBID}_Level_0.csv:ascii_csv

# 16 nodes
# spindle --location=/var/tmp/mcgurn4 srun -n576 /usr/workspace/mcgurn4/ablateOpt/ablate \
#    --input $FILE \
#    -yaml::environment::title dm1-$SLURM_JOBID \
#    -yaml::timestepper::domain::faces [280,60,60] \
#    -log_view :/p/lustre2/ubchrest/ablateInputsslabBurner.pmma.3D.scaling/ConvEst_Refinement_${SLURM_JOBID}_Level_1.csv:ascii_csv

# 128 nodes
# spindle --location=/var/tmp/mcgurn4 srun -n4608 /usr/workspace/mcgurn4/ablateOpt/ablate \
#    --input $FILE \
#    -yaml::environment::title dm2-$SLURM_JOBID \
#    -yaml::timestepper::domain::faces [560,120,120] \
#    -log_view :/p/lustre2/ubchrest/ablateInputsslabBurner.pmma.3D.scaling/ConvEst_Refinement_${SLURM_JOBID}_Level_2.csv:ascii_csv

# 1024 nodes
# spindle --location=/var/tmp/mcgurn4 srun -n36864 /usr/workspace/mcgurn4/ablateOpt/ablate \
#    --input $FILE \
#    -yaml::environment::title dm3-$SLURM_JOBID \
#    -yaml::timestepper::domain::faces [1120,240,240] \
#    -log_view :/p/lustre2/ubchrest/ablateInputsslabBurner.pmma.3D.scaling/ConvEst_Refinement_${SLURM_JOBID}_Level_3.csv:ascii_csv


echo 'Done'
