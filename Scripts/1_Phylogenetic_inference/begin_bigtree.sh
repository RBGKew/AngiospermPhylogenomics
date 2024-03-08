#!/bin/bash
#SBATCH --job-name=Begin_BigTree
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#BATCH --partition=short


mkdir -p 00_logs 01_seqs 02_aligned 03_cleaned 04_genetree 05_treeshrink 06_speciestree

SCRIPT_PATH="~/scripts/bigtree"

START_DIR=$(pwd -LP)
dataset=$(pwd | awk -F/ '{print $(NF)}')
prefix=$(echo $dataset | awk '{print substr($0,0,6)}')

touch ../${dataset}_running_$iteration.to

if [ "$iteration" == "i1" ]
  then
    seqs_dir="01_seqs_i1"
    aligned_dir="02_aligned_i1" 
    cleaned_dir="03_cleaned_i1" 
    genetree_dir="04_genetree_i1"
  else
    seqs_dir="06_seqs_i2"
    aligned_dir="07_aligned_i2"
    cleaned_dir="08_cleaned_i2" 
    genetree_dir="09_genetree_i2"
    reverse_string="--reverse $START_DIR/seqs_to_reverse.txt"
    drop_string="--drop_list $START_DIR/outliers_list_i1_clean.txt --dropped_only"
    awk '{gsub("_i1","");print}' $START_DIR/outliers_list_i1.txt > $START_DIR/outliers_list_i1_clean.txt 
    sh $SCRIPT_PATH/concat_misdirection_files.sh 01_seqs_i1/fix $START_DIR/seqs_to_reverse.txt
fi
mkdir -p 00_logs $seqs_dir/fix $aligned_dir $cleaned_dir $genetree_dir
cd $seqs_dir

python3 $SCRIPT_PATH/sequence_handler.py \
-f $START_DIR/$seq_file \
-v --by_gene \
--min_length 50 \
$reverse_string \
$drop_string

if [ "$iteration" == "i1" ]
  then 
    find . -name "g*.fasta" -type f -exec awk -v x=$min_seqs 'NR==x{exit 1}' {} \; -exec rm -f {} \;
fi

ls g*.fasta | cut -d "." -f 1 > ../selected_genes_$iteration.txt
n_genes=$(wc -l ../selected_genes_$iteration.txt | awk '{print $1}')

cd ../00_logs

MAFFT_JOBID=$(sbatch --job-name=${prefix}_${iteration}_mafft --mem=$iqtree_mem --partition=$partition --array=1-${n_genes}%$task_throttle $SCRIPT_PATH/ordinal_mafft_iqtree.sh $iteration ../selected_genes_$iteration.txt) \
&& cd $START_DIR && \
CHECH_MAFFT_JOBID=$(sbatch --mem=100M --job-name=${prefix}_${iteration}_check_mafft --dependency=afterany:${MAFFT_JOBID##* } --chdir=$START_DIR $SCRIPT_PATH/check_mafft_iqtree.sh $iteration ${MAFFT_JOBID##* } $astral_mem $partition ) && \
echo $MAFFT_JOBID > $START_DIR/${prefix}_mafft_${iteration}_jobid.txt && \
echo $CHECH_MAFFT_JOBID > $START_DIR/${prefix}_check_mafft_${iteration}_jobid.txt
