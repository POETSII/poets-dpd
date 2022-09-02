#!/bin/sh
#SBATCH --partition=amd
#SBATCH --mem=8000
#SBATCH --job-name=TEST_JOB
#SBATCH --ntasks-per-node=64

export LD_LIBRARY_PATH=/home/dbt1c21/packages/tbb2019_20180718oss/lib/intel64/gcc4.7

echo "Hello"

bin/benchmark_engine NaiveDPDEngineHalfStepTBB uniform-16

