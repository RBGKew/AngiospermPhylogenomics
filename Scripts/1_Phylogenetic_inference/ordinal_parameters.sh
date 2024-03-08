#! /bin/bash

### 

IQTREE_TASK_THROTTLE=30
IQTREE_MEM="200M"
IQTREE_PARTITION="medium"
TREESHRINK_MEM="1G"
TREESHRINK_PARTITION="medium"
ASTRAL_MEM="1G"
ASTRAL_PARTITION="medium"

SCRIPT_PATH="~/scripts/bigtree"
START_DIR=$(pwd -LP)
DATASET=$(pwd | awk -F/ '{print $(NF)}')
PREFIX=$(echo $DATASET | awk '{print substr($0,0,6)}')

LOG="$START_DIR/Log_$PREFIX.txt"

### FILES
SEQ_FILE="$START_DIR/${DATASET}_files.txt"
GENES_TO_ALIGN_FILE="$START_DIR/genes_to_align_only_$iteration.txt"
GENES_TO_TREE_FILE="$START_DIR/genes_to_tree_$iteration.txt"
TREE_JOBID_PATH="$START_DIR/${PREFIX}_iqtree_${iteration}_jobid.txt"

### Setting directories
LOG_DIR="$START_DIR/00_logs"
if [ "$iteration" == "i1" ]
  then
    SEQS_DIR="$START_DIR/01_seqs_i1"
    ALIGNED_DIR="$START_DIR/02_aligned_i1" 
    CLEANED_DIR="$START_DIR/03_cleaned_i1" 
    GENETREE_DIR="$START_DIR/04_genetree_i1"
  else
    SEQS_DIR="$START_DIR/06_seqs_i2"
    ALIGNED_DIR="$START_DIR/07_aligned_i2" 
    CLEANED_DIR="$START_DIR/08_cleaned_i2" 
    GENETREE_DIR="$START_DIR/09_genetree_i2"
fi
TREESHRINK_DIR="$START_DIR/05_treeshrink_i1"
FINAL_FILES_DIR="$START_DIR/10_trees_aligns"
ASTRAL_DIR="$START_DIR/11_astral"

