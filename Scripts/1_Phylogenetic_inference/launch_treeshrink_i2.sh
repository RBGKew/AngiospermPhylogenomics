#! /bin/bash 
#SBATCH --mem=100M

echo The script check_mafft_iqtree.sh was called as: 
echo $0 $@

iteration=$1
JobID=$2
astral_mem=$3
partition=$4

SCRIPT_PATH="~/scripts/bigtree"
START_DIR=$(pwd -LP)
dataset=$(pwd | awk -F/ '{print $(NF)}')
prefix=$(echo $dataset | awk '{print substr($0,0,6)}')

sacct_output=${prefix}_sacct_mafft_iqtree_${iteration}_J_$JobID.txt

iqtree_mem=$(cut -d "|" -f 12 $sacct_output | head -n 1)

TREESHRINK_JOBID=$(sbatch --mem=$astral_mem --partition=$partition --job-name=${prefix}_${iteration}_treeshrink --chdir=$START_DIR $SCRIPT_PATH/ordinal_treeshrink.sh $iteration $astral_mem )
BEGIN_I2_JOBID=$(sbatch --job-name=${prefix}_begin_i2 --dependency=afterok:${TREESHRINK_JOBID##* } --chdir=$START_DIR $SCRIPT_PATH/begin_ordinal_iteration.sh i2 ${dataset}_files.txt $partition 10 $iqtree_mem $astral_mem )
echo $TREESHRINK_JOBID > $START_DIR/${prefix}_treeshrink_${iteration}_jobid.txt
echo $BEGIN_I2_JOBID > $START_DIR/${prefix}_begin_i2_jobid.txt
