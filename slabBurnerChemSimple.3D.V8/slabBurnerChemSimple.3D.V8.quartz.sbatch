#!/bin/tcsh
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 12
#SBATCH -J parSolve34
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb


##### Load Required modules
module load gcc/10.2.1

##### Launch parallel job using srun
srun -n432 /g/g15/mcgurn4/ablateOpt/ablate --input /usr/workspace/ubchrest/mcgurn4/ablateInputs/slabBurnerChemSimple.3D.V8/slabBurnerChemSimple.3D.V8.yaml  -yaml::environment::title slab_$SLURM_JOBID -yaml::flow::mesh::options::dm_refine 5
echo 'Done'