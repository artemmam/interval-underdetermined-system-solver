#!/bin/bash
#SBATCH --partition=intel-gpu
#SBATCH --nodelist=g5500-1
echo Parallel version with 32 x 32 grid and 92 procs
/usr/local/bin/docker-id agorchakov:bio mpirun -n 92 --mca btl_vader_single_copy_mechanism none python3 -W ignore passive_rehabilitation_system.py -Nu 32 -Nv 10 -v1="-20,10" -v2="-80,10" --plotting --parallel

