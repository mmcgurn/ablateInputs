#!/bin/tcsh
#### See https://hpc.llnl.gov/training/tutorials/slurm-and-moab#LC

##### These lines are for Slurm
#SBATCH -N 8
#SBATCH -J 3dSlab
#SBATCH -t 00:60:00
#SBATCH -p pdebug
#SBATCH --mail-type=ALL
#SBATCH -A sunyb


##### Load Required modules
module load chrest/quartz
module load gcc/10.2.1

##### Launch parallel job using srun
srun -n288 /g/g15/mcgurn4/ablateOpt/ablate --input /p/lustre1/mcgurn4/ablateInputs/slabBurnerChem.3D.V13/slabBurnerChem.3D.V13.yaml  -yaml::environment::title slab_$SLURM_JOBID -yaml::flow::mesh::options::dm_refine 1
echo 'Done'