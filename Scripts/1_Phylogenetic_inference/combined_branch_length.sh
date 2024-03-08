#!/bin/bash
#
#SBATCH --job-name=IQTree_const

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=25G
#SATCH --partition=main

source ../.././parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

dataset=combined_$(pwd | awk -F/ '{print $(NF)}')
align_dir="aligns/"

outgroups="s17275"

iqtree2 \
-s $align_dir \
--prefix $dataset \
-g top*constraint.nwk \
-t top*constraint.nwk \
-o $outgroups \
-m GTR+F+R \
--mem 25G \
-T ${SLURM_CPUS_PER_TASK} \
> ./log-iqtree.out \
2> ./log-iqtree.err

