#!/bin/bash
#
#SBATCH --job-name=MAFFT_IQTREE
#SBATCH --output=00_logs/slurm-iqtree_%A_%a.out
#SBATCH --error=00_logs/slurm-iqtree_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

iteration=$1
source ./parameters.sh

gene=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' $GENES_TO_TREE_FILE)
echo $gene

cd $GENETREE_DIR

iqtree2 \
-s $CLEANED_DIR/${gene}_${iteration}_cleaned.fasta \
-m GTR+G \
-B 1000 \
-T auto \
--keep-ident \
-pre ${gene}_${iteration} \
2>> $LOG_DIR/log_${gene}_$iteration-iqtree.err
