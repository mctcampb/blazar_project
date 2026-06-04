#!/bin/bash
export PBS_O_WORKDIR=/data/euclid/scratch/campb951/blazar_project
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

module load python/3.12 
source ~/.env/marti_env/bin/activate
which python

python compute_wp.py 
