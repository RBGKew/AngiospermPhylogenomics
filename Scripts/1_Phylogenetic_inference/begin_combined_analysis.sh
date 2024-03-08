#!/bin/bash
#SBATCH --job-name=Combined_analysis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --mem=500M

source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

mkdir -p $LOG_DIR

# Copy backbone gene trees
if [ ! -d ${BACKBONE_TREES_DIR} ]; then
  mkdir -p ${BACKBONE_TREES_DIR}/1_input ${BACKBONE_TREES_DIR}/2_collapsed
  cp ${BACKBONE_ANALYSIS_DIR}/10_trees_aligns/*.treefile ${BACKBONE_TREES_DIR}/1_input
  for file in ${BACKBONE_TREES_DIR}/1_input/*.treefile; do
    gene=$(basename $file | cut -d "_" -f 1)
    $APPS/newick-utils-1.6/src/nw_ed $file 'i & b<80' o \
    > ${BACKBONE_TREES_DIR}/2_collapsed/${gene}_collapsed.tre
  done
  remove_outliers=$(sbatch --partition=short --mem=4G --wrap="Rscript ${SCRIPT_PATH}/remove_backbone_outliers.R \
    ${SCRIPT_PATH}/tips.csv \
    ${BACKBONE_TREES_DIR}/2_collapsed/ \
    _collapsed.tre \
    ${BACKBONE_TREES_DIR} \
    $LOG_DIR \
    ${BACKBONE_TREES_DIR} ")
  remove_outliers_jobID=${remove_outliers##* }
else
  echo "Backbone directory found. Skipping this step"
fi

# Congregate alignments
if [ ! -d $ORDINAL_ALIGNED_DIR ]; then
  mkdir -p $ORDINAL_ALIGNED_DIR
  cat $GENE_LIST | while read gene; do
    gene_destination="${ORDINAL_ALIGNED_DIR}/${gene}/"
    mkdir -p ${gene_destination}
    cp -v ${ORDINAL_ANALYSES_DIR}/*/10_trees_aligns/${gene}_*.fasta $gene_destination
    for file in ${gene_destination}/${gene}*.fasta; do
      echo >> $file
    done
  done
else
  echo "Ordinal aligments directory found. Skipping this step"
fi

# Congregate sequences not aligned
if [ ! -d ${ORDINAL_NOT_ALIGNED_DIR} ]; then
  mkdir -p ${ORDINAL_NOT_ALIGNED_DIR}/tmp/
  cp ${ORDINAL_ANALYSES_DIR}/*/10_trees_aligns/*_seqs_not_aligned.fasta ${ORDINAL_NOT_ALIGNED_DIR}/tmp/.
  find ${ORDINAL_NOT_ALIGNED_DIR}/tmp/ > ${ORDINAL_NOT_ALIGNED_DIR}/ordinal_not_aligned_files.txt 
  cp ${ORDINAL_ANALYSES_DIR}/Orphan/Orphan_files.txt ${ORDINAL_NOT_ALIGNED_DIR}/.
  cat ${ORDINAL_NOT_ALIGNED_DIR}/ordinal_not_aligned_files.txt \
    ${ORDINAL_NOT_ALIGNED_DIR}/Orphan_files.txt \
    > ${ORDINAL_NOT_ALIGNED_DIR}/all_files_seqs_not_aligned.txt
  mkdir -p ${ORDINAL_NOT_ALIGNED_SEQS_DIR}/fix 
  python3 $SCRIPT_PATH/sequence_handler.py \
    -f ${ORDINAL_NOT_ALIGNED_DIR}/all_files_seqs_not_aligned.txt \
    -d $ORDINAL_NOT_ALIGNED_SEQS_DIR \
    --by_gene

  mkdir -p ${ORDINAL_NOT_ALIGNED_ALIGNED_DIR}
  cat $GENE_LIST \
    | parallel $SCRIPT_PATH/combined_mafft_orphans.sh \
    {} \
    $ORDINAL_NOT_ALIGNED_SEQS_DIR \
    $ORDINAL_NOT_ALIGNED_ALIGNED_DIR \
    $LOG_DIR 

  mkdir -p ${ORDINAL_NOT_ALIGNED_REORIENTED_DIR}
  python3 $SCRIPT_PATH/sequence_handler.py \
    -f ${ORDINAL_NOT_ALIGNED_ALIGNED_DIR} \
    -d ${ORDINAL_NOT_ALIGNED_REORIENTED_DIR} \
    --in_type gene \
    --by_gene \
    --no_gap \
    -v
else
  echo "Unaligned sequences directory found. Skipping this step"
fi

if [ ! -d ${PRE_MERGE_DIR} ]; then
  echo "Launching arrays"
  mkdir -p ${PRE_MERGE_DIR} ${ALIGNED_DIR} ${CLEANED_ALIGNED_DIR}/tmp 
  sh ${SCRIPT_PATH}/combined_launch_arrays.sh $MAFFT_MEM
else
  echo "Checking arrays"
  sh ${SCRIPT_PATH}/combined_check_arrays.sh
fi

