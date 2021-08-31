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
srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_1_$SLURM_JOBID  -yaml::timestepper::name dmRefine1 -yaml::flow::mesh::options::dm_refine 1 -log_view :fv.$SLURM_JOBID.1.csv:ascii_csv

srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_2_$SLURM_JOBID  -yaml::timestepper::name dmRefine2 -yaml::flow::mesh::options::dm_refine 2 -log_view :fv.$SLURM_JOBID.2.csv:ascii_csv

srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_3_$SLURM_JOBID  -yaml::timestepper::name dmRefine3 -yaml::flow::mesh::options::dm_refine 3 -log_view :fv.$SLURM_JOBID.3.csv:ascii_csv

srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_4_$SLURM_JOBID  -yaml::timestepper::name dmRefine4 -yaml::flow::mesh::options::dm_refine 4 -log_view :fv.$SLURM_JOBID.4.csv:ascii_csv

srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_5_$SLURM_JOBID  -yaml::timestepper::name dmRefine5 -yaml::flow::mesh::options::dm_refine 5 -log_view :fv.$SLURM_JOBID.5.csv:ascii_csv

srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_6_$SLURM_JOBID  -yaml::timestepper::name dmRefine6 -yaml::flow::mesh::options::dm_refine 6 -log_view :fv.$SLURM_JOBID.6.csv:ascii_csv

srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_7_$SLURM_JOBID  -yaml::timestepper::name dmRefine7 -yaml::flow::mesh::options::dm_refine 7 -log_view :fv.$SLURM_JOBID.7.csv:ascii_csv

srun -n10800 ablate --input ../scalability.3D.V2.yaml  -yaml::environment::title slab_8_$SLURM_JOBID  -yaml::timestepper::name dmRefine8 -yaml::flow::mesh::options::dm_refine 8 -log_view :fv.$SLURM_JOBID.8.csv:ascii_csv
