#!/bin/sh

#SBATCH --partition=amd
#SBATCH --ntasks-per-node=64     # Tasks per node
#SBATCH --nodes=1                # Number of nodes requested

module load openmpi/4.1.1/amd-intel
module load lammps/2020/amd-intel

which lmp

bash run_dpd_water_mpi_helper.sh mpi_amd_n1_t64
