#!/bin/bash
#
#SBATCH --job-name=Ord_TreeShrink
#SBATCH --output=00_logs/slurm-treeShrink_%j.out
#SBATCH --error=00_logs/slurm-treeShrink_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --mem=100M

#echo The script ordinal_treeshrink.sh was called as: 
#echo $0 $@

iteration=$1

source ./parameters.sh

echo "" >> $LOG
date >> $LOG
echo "Beginning TreeShrink" >> $LOG

mkdir -p $TREESHRINK_DIR
cd $TREESHRINK_DIR

cat $GENETREE_DIR/*.treefile > ${PREFIX}_$iteration.trees

python $APPS/TreeShrink/run_treeshrink.py \
-t ${PREFIX}_$iteration.trees \
--centroid \
-m 'all-genes'

basename -a $GENETREE_DIR/*.treefile | cut -d "." -f 1 | awk -v var=$iteration '{gsub("_"var,""); print $0}' > $START_DIR/complete_gene_trees_$iteration.txt
paste -d ":" $START_DIR/complete_gene_trees_$iteration.txt ${PREFIX}_${iteration}_treeshrink/output.txt > $START_DIR/outliers_list_$iteration.txt

rm ../../${DATASET}_running_$iteration.to
touch ../../${DATASET}_${iteration}_done.to

echo "TreeShrink complete" >> $LOG
