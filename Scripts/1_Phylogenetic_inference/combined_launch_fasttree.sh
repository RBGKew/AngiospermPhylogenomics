#! /bin/bash
source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

array_list=""

squeue --me > tmp_squeue.txt


for i in $(seq 353); do
  gene=$(awk -v lineid=$i 'NR==lineid{print;exit}' $GENE_LIST)
  tree_fn="${GENE_TREES_DIR}/${gene}-fasttree.tre"
  if [ -s "$tree_fn" ]; then
    continue
  else
    running=$(grep "_$i " tmp_squeue.txt | wc -l)
    if [ "$running" -eq 0 ]; then
      array_list="${array_list}${i},"
    fi
  fi
done

array_list=${array_list::-1}
echo $array_list

echo sbatch \
 --mem=2G \
 --partition=medium \
  --array=$array_list%${FASTTREE_TASK_THROTTLE} \
  ${SCRIPT_PATH}/combined_fasttree.sh

rm tmp_squeue.txt
