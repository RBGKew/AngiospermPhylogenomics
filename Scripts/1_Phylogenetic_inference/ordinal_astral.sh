#!/bin/bash
#
#SBATCH --job-name=Ord_Astral
#SBATCH --output=00_logs/slurm-astral_%j.out
#SBATCH --error=00_logs/slurm-astral_%j.out
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium

source ./parameters.sh

echo ""
date >> $LOG
echo "Beginnig Astral" >> $LOG

mkdir -p $FINAL_FILES_DIR
mkdir -p $ASTRAL_DIR

sh $SCRIPT_PATH/colate_ordinal_files.sh -v

cd $ASTRAL_DIR
cat $FINAL_FILES_DIR/*.treefile > ${PREFIX}.trees

$APPS/newick-utils-1.6/src/nw_ed ${PREFIX}.trees 'i & b<=30' o > ${PREFIX}-BS30.trees

java -Xmx${ASTRAL_MEM} -jar $APPS/Astral/astral.5.7.3.jar \
-i ${PREFIX}-BS30.trees \
-o ${PREFIX}-BS30_sp.tre \
-t 2 \
2> ../00_logs/astral-BS30.log

cp ${PREFIX}-BS30_sp.tre ~/bigtree_results/v3/${DATASET}_astral_$(date +"%F").tre
rm ../../${DATASET}_running_i*.to
touch ../../${DATASET}_COMPLETED.to

echo "Astral analysis complete" >> $LOG

