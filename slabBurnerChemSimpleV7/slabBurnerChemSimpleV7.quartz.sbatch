#!/bin/tcsh
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 4
#SBATCH -J parSolve34
#SBATCH -t 00:60:00
#SBATCH -p pdebug
#SBATCH --mail-type=ALL
#SBATCH -A sunyb


##### Load Required modules
module load gcc/10.2.1

##### Launch parallel job using srun
srun -n144 /g/g15/mcgurn4/ablateOpt/ablate --input /usr/workspace/ubchrest/mcgurn4/ablateInputs/slabBurnerChemSimpleV7/slabBurnerChemSimpleV7.yaml  -yaml::environment::title slab_$SLURM_JOBID -yaml::flow::mesh::options::dm_refine 6
echo 'Done'