#!/bin/sh
#SBATCH --partition=amd
#SBATCH --mem=16000
#SBATCH --job-name=TEST_JOB
#SBATCH --ntasks-per-node=64

export LD_LIBRARY_PATH=/home/dbt1c21/packages/oneTBB-2019/build/linux_intel64_gcc_cc11.1.0_libc2.17_kernel3.10.0_release

module load gcc/11.1.0

echo "Hello"

bin/benchmark_engine avx2_dpd_engine  uniform-32

