#!/bin/sh

#SBATCH --ntasks-per-node=40     # Tasks per node
#SBATCH --nodes=1                # Number of nodes requested

module load lammps/2020/intel

mpiexec -np $SLURM_NTASKS lmp -in dpd_water_100x100x100_t1000.txt
