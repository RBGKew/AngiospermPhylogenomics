#!/bin/bash
#BATCH --job-name=Sequence_handler
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#BATCH --partition=short

echo The summary.sh was called as: 
echo $0 $@
SCRIPT_PATH="~/scripts/bigtree"

START_DIR=$(pwd -LP)
dataset=$(pwd | awk -F/ '{print $(NF)}')
prefix=$(echo $dataset | awk '{print substr($0,0,6)}')

python3 $SCRIPT_PATH/sequence_handler.py \
-f $START_DIR \
-v \
--in_type gene \
--summary \
--length \
--no_gap \
-o no_gap
