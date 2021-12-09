#!/bin/bash
#SBATCH --partition=intel-gpu
#SBATCH --nodelist=g5500-1
##echo Parallel version with $N x $N x $N grid and $P procs
echo Parallel version with 15 x 15 x 15 grid and 92 procs
/usr/local/bin/docker-id agorchakov:bio mpirun -n 92 --mca btl_vader_single_copy_mechanism none python3 -W ignore passive_regabilitation_system_3d.py -Nu 15 -Nv 5 -v1 15,60 -v2 45,150 -v3 15,45 --parallel --save_boxes


