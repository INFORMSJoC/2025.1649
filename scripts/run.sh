#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 1
export SLURM_EXACT=1
source /public1/soft/modules/module.sh
module load cmake/3.27.5
module load gcc/8.1.0
export PATH=/public1/home/t0s000936/Algorithm/build:$PATH
cd build
cmake -DCMAKE_CXX_FLAGS="-O3" -DCMAKE_CXX_COMPILER=/public1/soft/gcc/8.1.0/bin/g++ -DCMAKE_C_COMPILER=/public1/soft/gcc/8.1.0/bin/gcc ..
make
./CircleCovering C4 100 1 &
sleep 1
wait