#! /bin/bash

dataset=$(pwd | awk -F/ '{print $(NF)}')

mafft_mem=$(cat /mnt/shared/scratch/azuntini/private/bigtree/v1/ordinal/$dataset/*_sacct_mafft_iqtree_i1_J_*.txt | cut -d "|" -f 15 | sort -g | tail -n 1)
echo mafft MAX: $(echo "$mafft_mem / 2^20" | bc ) Mb

format=JobID,JobIDRaw,JobName,Partition,State,ExitCode,Start,End,Elapsed,ElapsedRaw,AllocCPUS,ReqMem,AveRSS,AveVMSize,MaxRSS,MaxVMSize
format_IFS=$(echo $format | awk 'gsub(","," ")')

ASTRAL_JOBID=$(cat /mnt/shared/scratch/azuntini/private/bigtree/v1/ordinal/$dataset/*_astral_i1_jobid.txt) 
echo $ASTRAL_JOBID

astral_mem=$(sacct -j ${ASTRAL_JOBID##* }  --noconvert --parsable --format=$format | grep "|batch|" | cut -d "|" -f 15 | sort -g | tail -n 1)
echo mafft MAX: $(echo "$astral_mem / 2^20" | bc ) Mb
