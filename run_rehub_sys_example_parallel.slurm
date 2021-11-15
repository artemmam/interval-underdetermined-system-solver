#!/bin/bash
#SBATCH --partition=intel-gpu
#SBATCH --nodelist=g5500-1
/usr/local/bin/docker-id agorchakov:bio mpirun -n 80 --mca btl_vader_single_copy_mechanism none python3 -W ignore passive_rehabilitation_system.py -Nu 80 -Nv 10 -v1_0 15 -v1_1 60 -v2_0 45 -v2_1 150

