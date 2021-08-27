#!/bin/tcsh
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 24
#SBATCH -J mile1
#SBATCH -t 4:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb


##### Load Required modules
module load chrest/release

##### Launch parallel job using srun
srun -n864 ablate --input /usr/workspace/ubchrest/mcgurn4/ablateInputs/scalability.2D.V0/dmRefine3/scalability.2D.V0.yaml  -yaml::environment::title slab_$SLURM_JOBID -yaml::flow::mesh::options::dm_refine 3 -log_view :fv.$SLURM_JOBID.csv:ascii_csv
echo 'Done'