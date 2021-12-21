#!/bin/bash
#SBATCH --partition=intel-gpu
#SBATCH --nodelist=g5500-1
##echo Parallel version with $N x $N x $N grid and $P procs
echo Parallel version with 15 x 15 x 15 grid and 92 procs
/usr/local/bin/docker-id agorchakov:bio mpirun -n 92 --mca btl_vader_single_copy_mechanism none python3 -W ignore passive_regabilitation_system_3d.py -Nu 52 -Nv 10 -v1="-20,10" -v2="-80,10" -v3="0,25" --parallel --save_boxes
