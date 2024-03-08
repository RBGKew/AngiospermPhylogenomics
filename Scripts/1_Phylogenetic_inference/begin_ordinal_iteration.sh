#!/bin/bash
#SBATCH --job-name=Begin_ordinal_iteration
#SBATCH --output=00_logs/slurm-begin_i_%j.out
#SBATCH --error=00_logs/slurm-begin_i_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#BATCH --partition=short

echo The begin_ordinal_iteration.sh was called as: 
echo $0 $@

iteration=$1

source ./parameters.sh

date >> $LOG
echo "The beginning $iteration iteration" >> $LOG
echo "" >> $LOG

touch ../${DATASET}_running_$iteration.to

if [ "$iteration" == "i2" ]
 then
    reverse_string="--reverse $START_DIR/seqs_to_reverse.txt"
    drop_string="--drop_list $START_DIR/outliers_list_i1.txt --dropped_only"
    sh $SCRIPT_PATH/concat_misdirection_files.sh 01_seqs_i1/fix $START_DIR/seqs_to_reverse.txt
    rm ../${DATASET}_i1_done.to
fi
mkdir -p $LOG_DIR $SEQS_DIR/fix $ALIGNED_DIR 
cd $SEQS_DIR

python3 $SCRIPT_PATH/sequence_handler.py \
-f $SEQ_FILE \
-v --by_gene \
--min_length 50 \
$reverse_string \
$drop_string

for file in *.fasta
  do gene=${file/.fasta}
  len=$(grep ">" $file | wc -l)
  if [ "$len" -ge "20" ]; then
    echo $gene >> $GENES_TO_TREE_FILE
  elif [ "$len" -ge "3" ] ; then
    echo $gene >> $GENES_TO_ALIGN_FILE
  else
    echo $gene >> ../genes_not_to_align_$iteration.txt
  fi
done

cd $START_DIR 

if [ -f $GENES_TO_ALIGN_FILE ]; then
  n_align=$(wc -l $GENES_TO_ALIGN_FILE | awk '{print $1}' )
  echo "$n_align gene(s) to align only"
  align_output="$LOG_DIR/slurm-mafft_${iteration}_%j.out"
  align=$(sbatch --output=$align_output --error=$align_output -J${PREFIX}_mafft_only --cpus-per-task=4 $SCRIPT_PATH/ordinal_mafft_only.sh $iteration )
  align_jobID=${align##* }
  align_dependency="--dependency=afterany:"$align_jobID
  echo "$n_align gene(s) to align only - JobID: $align_jobID " >> $LOG 
fi

if [ -f $GENES_TO_TREE_FILE ]; then
  n_genes=$(wc -l $GENES_TO_TREE_FILE | awk '{print $1}')
  echo "$n_genes gene(s) to tree"
  tree_output="$LOG_DIR/slurm-iqtree_${iteration}_%A_%a.out"
  check_tree_output="$LOG_DIR/slurm-check_${iteration}_%j.out"
  mkdir -p $CLEANED_DIR/tmp $GENETREE_DIR
  tree=$(sbatch --output=$tree_output --error=$tree_output --job-name=${PREFIX}_${iteration}_mafft_iqtree --mem=$IQTREE_MEM --partition=$IQTREE_PARTITION --array=1-${n_genes}%$IQTREE_TASK_THROTTLE $SCRIPT_PATH/ordinal_mafft_iqtree.sh $iteration ) && \
  tree_jobID=${tree##* } && \
  echo $tree_jobID > $TREE_JOBID_PATH && \
  check_tree=$(sbatch --output=check_tree_output --error=check_tree_output --mem=100M --job-name=${PREFIX}_${iteration}_check_mafft --dependency=afterany:$tree_jobID $SCRIPT_PATH/ordinal_check_mafft_iqtree.sh $iteration ) # $tree_JobID_path   $astral_mem $partition ) && \
  check_tree_jobID=${check_tree##* } &&\
  echo $check_tree_jobID  > ${PREFIX}_check_mafft_${iteration}_jobid.txt
  echo "$n_genes gene(s) to align and tree - JobID: $tree_jobID " >> $LOG 
else
  echo "No genes to tree estimation. Proceeding to ordinal files colation." >> $LOG 
  sbatch $align_dependency --mem=100M $SCRIPT_PATH/colate_ordinal_files.sh
  touch ../${DATASET}_COMPLETED.to
  rm ../${DATASET}_running_$iteration.to
fi
