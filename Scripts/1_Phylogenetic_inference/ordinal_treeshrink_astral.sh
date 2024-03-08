#!/bin/bash
#
#SBATCH --job-name=Ord_TreeShrink
#SBATCH --output=00_logs/slurm-treeShrink_%j.out
#SBATCH --error=00_logs/slurm-treeShrink_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --mem=4G

#echo The script ordinal_treeshrink.sh was called as: 
#echo $0 $@

iteration=$1

source ./parameters.sh

TREESHRINK_DIR=$START_DIR/00_TS
GENETREE_DIR=$START_DIR/04_genetree_i1

mkdir -p $TREESHRINK_DIR
cd $TREESHRINK_DIR

cat ~/scripts/bigtree/gene_list.txt |
  while read gene; do
    if [ -f $GENETREE_DIR/${gene}_i1.contree ]; then 
      cat $GENETREE_DIR/${gene}_i1.treefile >> ${PREFIX}_$iteration.pre.trees
    fi
  done

$APPS/newick-utils-1.6/src/nw_ed ${PREFIX}_$iteration.pre.trees 'i & b<=30' o > ${PREFIX}_$iteration.trees 

python $APPS/TreeShrink/run_treeshrink.py \
-t ${PREFIX}_$iteration.trees \
--centroid \
-m 'all-genes'

java -Xmx${ASTRAL_MEM} -jar $APPS/Astral/astral.5.7.3.jar \
-i backbone_i1_treeshrink/output.trees \
-o ${PREFIX}-BS30_sp.tre \
-t 2 
2> ../00_logs/astral-BS30.log
