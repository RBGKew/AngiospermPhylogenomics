#!/bin/bash
#
#SBATCH --job-name=MAFFT_parallel
#SBATCH --output=00_logs/slurm-mafft_%j.txt
#SBATCH --error=00_logs/slurm-mafft_%j.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M

iteration=$1
gene=$2

source ./paramaters.sh

cd $ALIGNED_DIR

mafft \
--retree 2 \
--maxiterate 1000 \
--adjustdirection \
--thread 1 \
$SEQS_DIR/${gene}.fasta \
> ./${gene}_${iteration}_aligned.fasta \
2> $LOG_DIR/log_${gene}_$iteration-mafft.out

grep ">_R_" ./${gene}_${iteration}_aligned.fasta > $SEQS_DIR/fix/${gene}_direction.txt
