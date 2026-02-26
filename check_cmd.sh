#!/bin/sh
#SBATCH --partition=amd
#SBATCH --mem=16000
#SBATCH --job-name=TEST_JOB
#SBATCH --ntasks-per-node=64
#SBATCH --nodes=1


module load gcc/11.1.0

echo "Hello"
echo "PATH=$PATH"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

lscpu

#ldd $(which povray)

#povray /home/dbt1c21/projects/poets-dpd/scratch/dpd_dmpci_grids/2022-06-21/2b6b/output/11/11.000010000.pov
