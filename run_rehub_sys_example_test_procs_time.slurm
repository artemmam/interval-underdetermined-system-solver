#!/bin/bash
#SBATCH --partition=intel-gpu
#SBATCH --nodelist=g5500-1
echo Parallel version with 32 x 32 grid and 1 proc
for i in {1..3}
do
    /usr/local/bin/docker-id agorchakov:bio python3 -W ignore passive_rehabilitation_system.py -Nu 32 -Nv 10 -v1_0 15 -v1_1 60 -v2_0 45 -v2_1 150 --record_time
done

