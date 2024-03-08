#!/bin/bash
#
#SBATCH --job-name=MAFFT_parallel
#
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M

gene=$1
ORDINAL_NOT_ALIGNED_SEQS_DIR=$2
ORDINAL_NOT_ALIGNED_ALIGNED_DIR=$3
LOG_DIR=$4

outfile="${ORDINAL_NOT_ALIGNED_ALIGNED_DIR}/${gene}.fasta"

mafft \
  --adjustdirection \
  --thread 1 \
  ${ORDINAL_NOT_ALIGNED_SEQS_DIR}/${gene}.fasta \
  > $outfile \
  2> ${LOG_DIR}/log_orphans_${gene}-mafft.out

grep ">_R_" $outfile > ${ORDINAL_NOT_ALIGNED_SEQS_DIR}/fix/${gene}_direction.txt
sed -i 's/>_R_/>/g' $outfile
