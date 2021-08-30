#!/bin/tcsh
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 300
#SBATCH -J mile1
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb

# make a temp dir to run ablate/petsc in
mkdir tmp_$SLURM_JOBID
cd tmp_$SLURM_JOBID

##### Load Required modules
module load  intel/19.0.4
module load chrest/release

##### Launch parallel job using srun
srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_$SLURM_JOBID  -yaml::timestepper::name dmRefine0 -yaml::flow::mesh::options::dm_refine 0 -log_view :fv.$SLURM_JOBID.csv:ascii_csv
echo 'Done'