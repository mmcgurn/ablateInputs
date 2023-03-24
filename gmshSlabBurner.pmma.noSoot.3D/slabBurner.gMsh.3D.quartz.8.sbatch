#!/bin/bash
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 8
#SBATCH -J gMshNoSoot
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

export DM_REFINE=0
export TITLE=lowG-gMsh-8n-dm$DM_REFINE-pmma-$SLURM_JOBID
export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.noSoot.3D/slabBurner3D.noSoot.lowG.3_8_23.yaml

# export DM_REFINE=0
# export TITLE=highG-gMsh-8n-dm$DM_REFINE-pmma-$SLURM_JOBID
# export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.noSoot.3D/slabBurner3D.noSoot.highG.3_8_23.yaml

# export DM_REFINE=1
# export TITLE=lowG-gMsh-8n-dm$DM_REFINE-pmma-$SLURM_JOBID
# export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.noSoot.3D/slabBurner3D.noSoot.lowG.3_8_23.yaml

# export DM_REFINE=1
# export TITLE=highG-gMsh-8n-dm$DM_REFINE-pmma-$SLURM_JOBID
# export FILE=/p/lustre2/ubchrest/ablateInputs/gmshSlabBurner.pmma.noSoot.3D/slabBurner3D.noSoot.highG.3_8_23.yaml

##### Launch parallel job using srun
spindle --location=/var/tmp/mcgurn4 srun -n288 /usr/workspace/mcgurn4/ablateOpt/ablate \
   --input $FILE \
   -yaml::environment::title $TITLE \
   -yaml::timestepper::domain::options::dm_refine $DM_REFINE

echo 'Done'
