#! /bin/bash

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

