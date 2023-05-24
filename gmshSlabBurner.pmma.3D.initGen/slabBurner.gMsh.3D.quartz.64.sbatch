#!/bin/bash
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 64
#SBATCH -J gMsh
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb
#SBATCH --mail-user=mtmcgurn@buffalo.edu

##### Load Required modules
# gcc
module load clang/14.0.6-magic
module load cmake/3.25.2

# Load PETSC ENV
export PETSC_DIR="/usr/workspace/mcgurn4/petsc"
export PETSC_ARCH="arch-ablate-opt" # arch-ablate-debug or arch-ablate-opt
export PKG_CONFIG_PATH="${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$PKG_CONFIG_PATH"
export HDF5_ROOT="${PETSC_DIR}/${PETSC_ARCH}"  
# Include the bin directory to access mpi commands
export PATH="${PETSC_DIR}/${PETSC_ARCH}/bin:$PATH"

# Make a temp directory so that tchem has a place to vomit its files
mkdir tmp_$SLURM_JOBID
cd tmp_$SLURM_JOBID

export DM_REFINE=0
export TITLE=lowG-gMsh-64n-dm$DM_REFINE-HF-$SLURM_JOBID
export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.3D.initGen/slabBurner3D.lowG.3_8_23.HF.yaml


# export DM_REFINE=0
# export TITLE=lowG-gMsh-64n-dm$DM_REFINE-ls-$SLURM_JOBID
# export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.3D.initGen/slabBurner3D.lowG.3_8_23.LS.yaml

# export DM_REFINE=0
# export TITLE=lowG-gMsh-64n-dm$DM_REFINE-pmma-$SLURM_JOBID
# export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.3D.initGen/slabBurner3D.lowG.3_8_23.yaml



##### Launch parallel job using srun
srun -n2304 /usr/workspace/mcgurn4/ablateOpt/ablate \
   --input $FILE \
   -yaml::environment::title $TITLE \
   -yaml::timestepper::domain::options::dm_refine $DM_REFINE -build_twosided redscatter

echo 'Done'
