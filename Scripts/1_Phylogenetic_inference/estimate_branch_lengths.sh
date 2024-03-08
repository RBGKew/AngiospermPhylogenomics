#!/bin/bash
#
#SBATCH --job-name=BrLens
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium


n_top_genes=$1

source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

mkdir -p $CALIBRATION_DIR
cd CALIBRATION_DIR

${SCRIPT_PATH}/./prepare_branch_lenght_alignments.R \
    -f ${SORTADATE_DIR}/combined_gg.txt \
    -g ${n_top_genes} \
    -s ${SORTADATE_DIR}/astral.tre \
    -t ".tre" \
    -a ${SORTADATE_DIR}/aligns \
    -b ".fasta"

