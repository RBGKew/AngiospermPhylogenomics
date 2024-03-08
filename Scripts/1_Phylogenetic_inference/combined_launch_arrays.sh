#! /bin/bash

source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

array_mem=$1

if [ -z "$2" ]; then
  array=1-353
else
  array=$2
fi

if [ -d "${GENE_TREES_DIR}" ]; then
  mode="fasttree"
  check_output=${FASTTREE_JOBID_PATH}
else
  mode="mafft"
  check_output=${MAFFT_JOBID_PATH}
fi

if [ "$mode" == "mafft" ]; then
  array_job=$(sbatch \
  --mem=${array_mem} \
  --partition=$MAFFT_PARTITION \
  --array=${array}%${MAFFT_TASK_THROTTLE} \
  ${SCRIPT_PATH}/combined_mafft.sh)
else
  array_job=$(sbatch \
    --mem=${array_mem} \
    --partition=$FASTTREE_PARTITION \
    --array=${array}%${FASTTREE_TASK_THROTTLE} \
    ${SCRIPT_PATH}/combined_fasttree.sh)
fi
array_jobID=${array_job##* }
echo ${array_jobID} >> ${FASTTREE_JOBID_PATH} && \
check_array=$(sbatch \
  --dependency=afterany:${array_jobID} \
  ${SCRIPT_PATH}/combined_check_arrays.sh )
check_array_jobID=${check_array##* } &&\
echo ${check_array_jobID} >> check_${mode}_jobid.txt
