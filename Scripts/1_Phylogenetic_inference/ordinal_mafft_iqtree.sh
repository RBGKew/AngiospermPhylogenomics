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

cd $ALIGNED_DIR

wc $SEQS_DIR/${gene}.fasta

mafft \
--retree 2 \
--maxiterate 1000 \
--adjustdirection \
--thread 1 \
$SEQS_DIR/${gene}.fasta \
> ./${gene}_${iteration}_aligned.fasta \
2> $LOG_DIR/log_${gene}_$iteration-mafft.out

mkdir -p $CLEANED_DIR/tmp

cd $CLEANED_DIR

java -jar $APPS/phyutility/phyutility.jar \
-clean 0.10 \
-in $ALIGNED_DIR/${gene}_${iteration}_aligned.fasta \
-out ./tmp/${gene}_${iteration}_cleaned.fasta

grep ">_R_" ./tmp/${gene}_${iteration}_cleaned.fasta > $SEQS_DIR/fix/${gene}_direction.txt
sed -i 's/>_R_/>/g' ./tmp/${gene}_${iteration}_cleaned.fasta

sed '/^-*$/d' ./tmp/${gene}_${iteration}_cleaned.fasta | grep -v --no-group-separator -B1 '^>' > ${gene}_${iteration}_cleaned.fasta

cd $GENETREE_DIR

iqtree2 \
-s $CLEANED_DIR/${gene}_${iteration}_cleaned.fasta \
-m GTR+G \
-B 1000 \
-T auto \
--keep-ident \
--redo \
-pre ${gene}_${iteration} \
2> $LOG_DIR/log_${gene}_$iteration-iqtree.err
