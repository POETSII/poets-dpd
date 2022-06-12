#!/bin/sh

#SBATCH --ntasks-per-node=40     # Tasks per node
#SBATCH --nodes=1                # Number of nodes requested

module load lammps/2020/intel

bash run_dpd_water_mpi_helper.sh mpi_intel_n1_t40
