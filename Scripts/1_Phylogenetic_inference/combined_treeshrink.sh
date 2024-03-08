#!/bin/bash
#
#SBATCH --job-name=TreeShrink
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=18G
#SBATCH --partition=long
#SBATCH --output=00_logs/TreeShrink.%A.log

source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

mkdir -p $TREESHRINK_DIR

cat $GENE_LIST | while read gene; do
  align_path="${CLEANED_ALIGNED_DIR}/${gene}_cleaned.fasta"
  tree_path="${GENE_TREES_DIR}/${gene}-fasttree.tre"
  if [ -s "$align_path" ] && [ -s "$tree_path" ]; then
    dest_dir="${TREESHRINK_DIR}/genes/$gene"
    mkdir -p $dest_dir
    cp $align_path $dest_dir/input.fasta
    $APPS/newick-utils-1.6/src/nw_ed $tree_path 'i & b<=0.1' o > $dest_dir/input.tre
  else
    echo $gene >> "${TREESHRINK_DIR}/missing_genes.txt"
  fi
done

cd $TREESHRINK_DIR
pwd

python $APPS/TreeShrink/run_treeshrink.py \
  -i ./genes/  \
  -t input.tre \
  -a input.fasta \
  --centroid \
  -m 'all-genes'
