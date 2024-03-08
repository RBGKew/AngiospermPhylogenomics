#! /bin/bash

#SBATCH --output=00_logs/slurm-check_array_%j.out
#SBATCH --error=00_logs/slurm-check_array_%j.out
#SBATCH --mem=100M

source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

if [ -d "${FASTTREE_JOBID_PATH}" ]; then
  mode="fasttree"
  jobID=$(tail -n 1 $FASTTREE_JOBID_PATH)
else
  mode="mafft"
  jobID=$(tail -n 1 $MAFFT_JOBID_PATH)
fi

echo "" >> $LOG
date >> $LOG
echo "Checking $mode" >> $LOG

sacct_format="JobID,JobIDRaw,JobName,Partition,State,ExitCode,Start,End,Elapsed,ElapsedRaw,AllocCPUS,ReqMem,AveRSS,AveVMSize,MaxRSS,MaxVMSize"
format_IFS=$(echo $sacct_format | awk 'gsub(","," ")')

sacct_output=${START_DIR}/sacct_${mode}_J_$jobID.txt
sacct -j $jobID --noconvert --parsable --format=$sacct_format | grep "|batch|" > $sacct_output

grep -v "|COMPLETED|" $sacct_output > sacct_batch.tmp

if [ -s sacct_batch.tmp ]; then
# RE DO  TASKS THAT FAILED IN THE PREVIOUS ARRAY
  highest_mem=$(cat $sacct_output | cut -d "|" -f 15 | sort -r | head -n 1)
  new_mem=$(echo "2 * $highest_mem / 2^20" | bc )
  echo $new_mem
  if [ "$new_mem" -lt "400" ]; then
    new_mem=2000
  fi

  echo "Some samples failed"
  cat sacct_batch.tmp | \
    while IFS="|" read -r $format_IFS ; do
      echo $JobID,$ReqMem,$MaxVMSize
      taskID=$(echo $JobID | cut -d "_" -f 2 | cut -d "." -f 1)
      echo $taskID >> tasks.tmp
      gene=$(awk "NR==$TaskID" $GENE_LIST)
    done
  tasks=$(paste -s -d, tasks.tmp)
  n_tasks=$(wc -l tasks.tmp | awk '{print $1'})

  sh ${SCRIPT_PATH}/combined_launch_arrays.sh $new_mem $tasks
  rm tasks.tmp

else
# ALL TASKS IN THE PREVIOUS ARRAY ARE COMPLETE
  if [ "$mode" == "mafft" ]; then
    mkdir -p ${GENE_TREES_DIR}/constraints
    sh ${SCRIPT_PATH}/combined_prepare_constraints.sh
    sh ${SCRIPT_PATH}/combined_launch_arrays.sh $FASTTREE_MEM
  else
    treeshrink=$(sbatch \
      --mem=$TREESHRINK_MEM \
      --partition=$TREESHRINK_PARTITION \
      --job-name=${PREFIX}_${iteration}_treeshrink \
      $SCRIPT_PATH/ordinal_treeshrink.sh $iteration ) && \
    treeshrink_jobID=${treeshrink##* } && \
    echo $treeshrink_jobID >> $START_DIR/treeshrink_jobid.txt

    ulimit -c unlimited
    ulimit -s unlimited
    ulimit -H unlimited
    astral=$(sbatch \
      --mem=$ASTRAL_MEM \
      --partition=$ASTRAL_PARTITION \
      --dependency=afterok:${treeshrink_jobID} \
      ${SCRIPT_PATH}/combined_astral.sh )
    astral_jobID=${astral##* }
    echo $astral_jobID >> $START_DIR/astral_jobid.txt
    echo "Second iteration complete. Proceeding to Astral - JobID: $astral_jobID" >> $LOG 
  fi
fi

rm sacct_batch.tmp
