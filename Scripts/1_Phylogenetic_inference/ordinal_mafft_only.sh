#!/bin/bash
#
#SBATCH --job-name=MAFFT_parallel
#
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M

iteration=$1

source ./parameters.sh

SLURM_CPUS_PER_TASK=4

function mafft_only() {  
  iteration=$1
  gene=$2
  seqs_path=$3

  echo $gene
  mafft \
  --retree 2 \
  --maxiterate 1000 \
  --adjustdirection \
  --thread 1 \
  $seqs_path/${gene}.fasta \
  > ./${gene}_${iteration}_aligned.fasta \
  2> ../00_logs/log_${gene}_$iteration-mafft.out

  grep ">_R_" ./${gene}_${iteration}_aligned.fasta > $seqs_path/fix/${gene}_direction.txt
  sed -i 's/>_R_/>/g' ./${gene}_${iteration}_aligned.fasta
}

export -f mafft_only

cd $ALIGNED_DIR

cat $GENES_TO_ALIGN_FILE | parallel --jobs ${SLURM_CPUS_PER_TASK} mafft_only $iteration {} $SEQS_DIR
