#!/bin/bash
#SBATCH --partition=intel-gpu
#SBATCH --nodelist=g5500-1
for i in 8 16 32
do
    for j in 4 8 16 32 64   
    do
        echo Parallel version with $i x $i grid and $j procs
        /usr/local/bin/docker-id agorchakov:bio mpirun -n $j --mca btl_vader_single_copy_mechanism none python3 -W ignore passive_rehabilitation_system.py -Nu $i -Nv 10 -v1_0 15 -v1_1 60 -v2_0 45 -v2_1 150 --parallel
    done
    echo Sequential version with $i x $i grid
    /usr/local/bin/docker-id agorchakov:bio python3 -W ignore passive_rehabilitation_system.py -Nu $i -Nv 10 -v1_0 15 -v1_1 60 -v2_0 45 -v2_1 150
done

