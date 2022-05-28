#!/bin/sh
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1

module load lammps/2018/cuda

echo
echo  ==============================
echo

mpirun -np 16 lmp -sf gpu -pk gpu 4 -in dpd_water_100x100x100_t1000.txt

echo
echo  ==============================
echo

mpirun -np 16 lmp -sf gpu -pk gpu 2 -in dpd_water_100x100x100_t1000.txt

echo
echo  ==============================
echo

mpirun -np 16 lmp -sf gpu -pk gpu 1 -in dpd_water_100x100x100_t1000.txt

#mpiexec -np $SLURM_NTASKS lmp -in dpd_water_100x100x100_t1000.txt
