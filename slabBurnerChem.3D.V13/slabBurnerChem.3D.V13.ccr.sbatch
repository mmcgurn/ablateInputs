#!/bin/sh
#SBATCH --partition=general-compute --qos=supporters
#SBATCH --time=72:00:00
#SBATCH --nodes=40
#SBATCH --constraint=NIH
#SBATCH --ntasks-per-node=40
#SBATCH --mem=192000
#SBATCH --job-name="slb2D"
#SBATCH --mail-user=mtmcgurn@buffalo.edu
#SBATCH --mail-type=END

# Print the current environment
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

# Load the required modules
module load intel/19.5 hdf5/1.12.0-mpi mkl/2019.5 gcc/10.2.0 cmake/3.17.1
 
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
srun -n $SLURM_NPROCS  /projects/academic/chrest/mtmcgurn/ablateOpt/ablate --input /panasas/scratch/grp-ped3/mtmcgurn/ablateInputs/slabBurnerChem.3D.V13/slabBurnerChem.3D.V13.yaml  -yaml::environment::title slab_$SLURM_JOBID -yaml::flow::mesh::options::dm_refine 1
echo "End Time " `date +%s`
