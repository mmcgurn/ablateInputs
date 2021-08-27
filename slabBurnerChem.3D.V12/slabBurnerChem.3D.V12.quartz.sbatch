#!/bin/tcsh
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 24
#SBATCH -J 3dSlab
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb


##### Load Required modules
module load gcc/10.2.1

##### Launch parallel job using srun
srun -n864 /g/g15/mcgurn4/ablateOpt/ablate --input /usr/workspace/ubchrest/mcgurn4/ablateInputs/slabBurnerChem.3D.V12/slabBurnerChem.3D.V12.yaml  -yaml::environment::title slab_$SLURM_JOBID -yaml::flow::mesh::options::dm_refine 4
echo 'Done'