#!/bin/sh
#SBATCH --partition=ped3 --qos=ped3
#SBATCH --clusters=faculty
#SBATCH --time=48:00:00
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10000
#SBATCH --job-name="slabFlow2D"
#SBATCH --mail-user=mtmcgurn@buffalo.edu
#SBATCH --mail-type=END

# Print the current environment
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

# Load the required modules
module load hdf5/1.12.0-mpi 
module load gcc/10.2.0
module load cmake/3.17.1

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
srun -n $SLURM_NPROCS  /projects/academic/ped3/mtmcgurn/ablate/opt/ablate --input /panasas/scratch/grp-ped3/mtmcgurn/ablateInputs/slabBurner2D/slabBurner2D.yaml
echo "End Time " `date +%s`
