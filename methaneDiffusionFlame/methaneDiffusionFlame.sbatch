#!/bin/sh
#SBATCH --partition=ped3 --qos=ped3
#SBATCH --clusters=faculty
#SBATCH --time=48:00:00
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10000
#SBATCH --job-name="slbHot"
#SBATCH --mail-user=mtmcgurn@buffalo.edu
#SBATCH --mail-type=END

# Print the current environment
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

# Load the required modules
module load chrest/release
module load petsc-chrest/release

# The initial srun will trigger the SLURM prologue on the compute nodes.
NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo NPROCS=$NPROCS

# The PMI library is necessary for srun
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Tell the tests what mpi command to use
export TEST_MPI_COMMAND=srun

# change to your build directory
cd debug
echo "current directory ="$PWD

# Run all tests
echo "Start Time " `date +%s`
# srun -n $SLURM_NPROCS valgrind --tool=memcheck -q --num-callers=20 --suppressions=/projects/academic/ped3/petsc/lib/petsc/bin/maint/petsc-val.supp --log-file=valgrind.log.%p  /projects/academic/ped3/mtmcgurn/ablate/opt/ablate --input /panasas/scratch/grp-ped3/mtmcgurn/ablateInputs/slabBurnerChem2D/slabBurnerChem2D.yaml
srun -n $SLURM_NPROCS  /projects/academic/chrest/mtmcgurn/ablateOpt/ablate --input /panasas/scratch/grp-ped3/mtmcgurn/ablateInputs/methaneDiffusionFlame.yaml -yaml::timestepper::arguments::ts_adapt_safety .1 -yaml::environment::title methaneDiffusionFlame_$SLURM_JOBID

echo "End Time " `date +%s`
