#!/bin/sh
#SBATCH -N 1
#SBATCH -J 2dSlab
#SBATCH -t 24:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A sunyb
#SBATCH --mail-user=mtmcgurn@buffalo.edu


module purge
module load clang/14.0.6-magic
module load cmake/3.25.2


export PETSC_DIR="/usr/workspace/mcgurn4/petsc" #UPDATE to the real path of petsc
export PETSC_ARCH="arch-ablate-opt"
export PKG_CONFIG_PATH="${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$PKG_CONFIG_PATH"
export HDF5_ROOT="${PETSC_DIR}/${PETSC_ARCH}"  
export PATH="${PETSC_DIR}/${PETSC_ARCH}/bin:$PATH"


srun -n36 /usr/workspace/mcgurn4/ablateOpt/ablate --input /p/lustre2/mcgurn4/ablateInputs/ablateCoda/slabs/slabBurner.2D.1.yaml
