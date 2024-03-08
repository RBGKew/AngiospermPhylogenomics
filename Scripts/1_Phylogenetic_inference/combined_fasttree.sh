#!/bin/bash
#
#SBATCH --job-name=FastTree
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=long
#SBATCH --output=00_logs/slurm-Fasttree.%A_%a.log

source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

gene=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' $GENE_LIST)

cd ${GENE_TREES_DIR}

fasttree \
  -nt \
  -gtr \
  -gamma \
  -pseudo \
  -spr 6 \
  -mlacc 3 \
  -slownni \
  -constraints ${CONSTRAINTS_DIR}/${gene}_constraint.txt \
  < ${CLEANED_ALIGNED_DIR}/${gene}_cleaned.fasta \
  > ${gene}-fasttree.tre \
  2> ${gene}-fasttree.txt
