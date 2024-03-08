#! /bin/bash

source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

mkdir -p ${CONSTRAINTS_DIR}

Rscript \
  ${SCRIPT_PATH}/prune_tree_to_alignment.R \
  ${CLEANED_ALIGNED_DIR} \
  "_cleaned.fasta" \
  ${BACKBONE_TREES_DIR} \
  "_backbone.tre" \
  ${CONSTRAINTS_DIR} \
  "_pruned_backbone.tre"

cat $GENE_LIST | \
  while read gene; do
  constraint_fn="${CONSTRAINTS_DIR}/${gene}_constraint.txt"
  if [ -f "$constraint_fn" ]; then
    continue
  fi
    perl \
      ${SCRIPT_PATH}/TreeToConstraints.pl \
      < ${CONSTRAINTS_DIR}/${gene}_pruned_backbone.tre \
      > ${constraint_fn}
  done
