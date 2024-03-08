#! /bin/bash

#SBATCH --output=00_logs/slurm-check_iqtree__%j.out
#SBATCH --error=00_logs/slurm-check_iqtree_%j.out


# The jobid is submitted with
iteration=$1

source ./parameters.sh

echo "" >> $LOG
date >> $LOG
echo "Checking IQ-Tree" >> $LOG

jobID=$(tail -n 1 $TREE_JOBID_PATH)

sacct_format="JobID,JobIDRaw,JobName,Partition,State,ExitCode,Start,End,Elapsed,ElapsedRaw,AllocCPUS,ReqMem,AveRSS,AveVMSize,MaxRSS,MaxVMSize"
format_IFS=$(echo $sacct_format | awk 'gsub(","," ")')

sacct_output=$START_DIR/${PREFIX}_sacct_mafft_iqtree_${iteration}_J_$jobID.txt
sacct -j $jobID --noconvert --parsable --format=$sacct_format | grep "|batch|" > $sacct_output

grep -v "|COMPLETED|" $sacct_output > sacct_batch.tmp

if [ -s sacct_batch.tmp ] && [ "$(wc -l $TREE_JOBID_PATH | awk '{print $1}' )" -lt "4" ] ; then
  highest_mem=$(cat $sacct_output | cut -d "|" -f 15 | sort -r | head -n 1)
  new_mem=$(echo "1.5 * $highest_mem / 2^20" | bc )
  echo $new_mem
  if [ "$new_mem" -lt "400" ]; then
    new_mem=400
  fi

  echo "Some samples failed"
  cat sacct_batch.tmp | while IFS="|" read -r $format_IFS
    do echo $JobID,$ReqMem,$MaxVMSize
    taskID=$(echo $JobID | cut -d "_" -f 2 | cut -d "." -f 1)
    echo $taskID >> tasks.tmp
    gene=$(awk "NR==$TaskID" ../genes_to_tree_$iteration.txt)
  done
  tasks=$(paste -s -d, tasks.tmp)
  n_tasks=$(wc -l tasks.tmp | awk '{print $1'})
  tree_output="$LOG_DIR/slurm-tree_%A_%a.out" && \
  tree=$(sbatch --output=$tree_output --error=$tree_output --job-name=${PREFIX}_${iteration}_mafft_iqtree --mem=$new_mem --partition=$IQTREE_PARTITION --array=$tasks%$IQTREE_TASK_THROTTLE $SCRIPT_PATH/ordinal_mafft_iqtree.sh $iteration ) && \
  tree_jobID=${tree##* } && \
  echo $tree_jobID >> $TREE_JOBID_PATH && \
  check_tree_output="$LOG_DIR/slurm-check_%j.out" && \
  check_tree=$(sbatch --output=$check_tree_output --error=$check_tree_output --mem=100M --job-name=${PREFIX}_${iteration}_check_mafft --dependency=afterany:$tree_jobID $SCRIPT_PATH/ordinal_check_mafft_iqtree.sh $iteration ) # $tree_JobID_path ) #$astral_mem $partition ) && \
  check_tree_jobID=${check_tree##* } && \
  echo $check_tree_jobID  > $START_DIR/${PREFIX}_check_mafft_${iteration}_jobid.txt
  echo "$n_tasks gene(s) to redo mafft and IQTree  - JobID: $TREE_JOBID" >> $LOG 
  rm tasks.tmp
else
  if [ "$iteration" == "i1" ]; then
    iqtree_mem=$(cut -d "|" -f 12 $sacct_output | head -n 1)
    echo "All_samples_worked. Proceeding to TreeShrink"
    treeshrink_output="$LOG_DIR/slurm-treeshrink_${iteration}_%j.out"
    treeshrink=$(sbatch --output=$treeshrink_output --error=$treeshrink_output --mem=$TREESHRINK_MEM --partition=$TREESHRINK_PARTITION --job-name=${PREFIX}_${iteration}_treeshrink $SCRIPT_PATH/ordinal_treeshrink.sh $iteration ) #$astral_mem )
    treeshrink_jobID=${treeshrink##* }
    echo $treeshrink_jobID > $START_DIR/${PREFIX}_treeshrink_${iteration}_jobid.txt
    begin_i2_output="$LOG_DIR/slurm-begin_i2_${iteration}_%j.out"
    begin_i2=$(sbatch --output=$begin_i2_output --error=$begin_i2_output --job-name=${PREFIX}_begin_i2 --dependency=afterok:$treeshrink_jobID $SCRIPT_PATH/begin_ordinal_iteration.sh i2 ) #$partition 10 $iqtree_mem $astral_mem )
    begin_i2_jobID=${begin_i2##* }
    echo $begin_i2_jobID > $START_DIR/${PREFIX}_begin_i2_jobid.txt
    echo "Proceeding to TreeShrink - JobID: $treeshrink_jobID" >> $LOG 
    echo "Iteration 2 will continue on JobID: $begin_i2_jobID " >> $LOG 
  else
    echo "Second iteration complete. Proceeding to Astral"
    colate_output="$LOG_DIR/slurm-colate_%j.out"
    colate=$(sbatch $align_dependency --output=$colate_output --error=$colate_output --mem=100M $SCRIPT_PATH/colate_ordinal_files.sh -v)
    colate_jobID=${colate##* }
    astral_output="$LOG_DIR//slurm-astral_%j.out"
    astral=$(sbatch --output=$astral_output --error=$astral_output --mem=$ASTRAL_MEM --partition=$ASTRAL_PARTITION --job-name=${PREFIX}_${iteration}_astral --dependency=afterany:$colate_jobID $SCRIPT_PATH/ordinal_astral.sh )
    astral_jobID=${astral##* }
    echo $astral_jobID > $START_DIR/${PREFIX}_astral_${iteration}_jobid.txt
    echo "Second iteration complete. Proceeding to Astral - JobID: $astral_jobID" >> $LOG 
  fi
fi

rm sacct_batch.tmp
