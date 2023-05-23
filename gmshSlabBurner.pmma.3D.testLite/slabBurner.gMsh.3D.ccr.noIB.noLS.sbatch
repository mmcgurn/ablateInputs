#!/bin/sh
#SBATCH --account=chrest
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --job-name="sbLite"
#SBATCH --mail-user=mtmcgurn@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --requeue

# Print the current environment
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

# Load the required modules
module unload intel-mpi/2020.2
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

# export DM_REFINE=0
# export TITLE=lowG-gMsh-dm$DM_REFINE-baseline-$SLURM_JOBID
# export FILE=/panasas/scratch/grp-chrest/mtmcgurn/ablateInputs/gmshSlabBurner.pmma.3D.testLite/slabBurner3D.lowG.defaultLim.yaml
# 
export DM_REFINE=0
export TITLE=lowG-gMsh-dm$DM_REFINE-defaultLimt-$SLURM_JOBID
export FILE=/panasas/scratch/grp-chrest/mtmcgurn/ablateInputs/gmshSlabBurner.pmma.3D.testLite/slabBurner3D.lowG.cold.noLS.yaml
# 
# 
# export DM_REFINE=0
# export TITLE=lowG-gMsh-dm$DM_REFINE-cold-$SLURM_JOBID
# export FILE=/panasas/scratch/grp-chrest/mtmcgurn/ablateInputs/gmshSlabBurner.pmma.3D.testLite/slabBurner3D.lowG.cold.yaml

# Run all tests
srun -n $SLURM_NPROCS  /projects/academic/chrest/mtmcgurn/ablateOpt/ablate \
  --input $FILE \
-yaml::environment::title $TITLE \
  -yaml::timestepper::domain::options::dm_refine $DM_REFINE


 