#! /bin/bash

if [ ! -e parameters.sh ]; then
  cp ~/scripts/bigtree/ordinal_parameters.sh ./parameters.sh
fi

source ./parameters.sh

mkdir -p $LOG_DIR

begin_i1_output="$LOG_DIR/slurm-begin_i1_%j.out"
sbatch -o$begin_i1_output -e$begin_i1_output $SCRIPT_PATH/begin_ordinal_iteration.sh i1
