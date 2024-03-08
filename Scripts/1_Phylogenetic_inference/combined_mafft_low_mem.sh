#!/bin/bash
#
#SBATCH --job-name=MAFFT_parallel
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=13000M
#SBATCH --partition=long
#SBATCH --output=00_logs/MAFFT_merge.%A_%a.log


source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh

gene=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' $GENE_LIST_SORTED)

cd ${PRE_MERGE_DIR}

not_aligned_seqs=${ORDINAL_NOT_ALIGNED_REORIENTED_DIR}/${gene}.fasta

if [ ! -f "$not_aligned_seqs" ]; then
  not_aligned_seqs=""
fi

gene_aligns="${ORDINAL_ALIGNED_DIR}/${gene}/${gene}_*.fasta"

cat $gene_aligns $not_aligned_seqs > ${gene}_pre_align.fasta

sh $SCRIPT_PATH/makemergetable.sh $gene_aligns > ${gene}_subMSAtable.txt

mafft \
  --maxiterate 10 \
  --merge ${gene}_subMSAtable.txt \
  --thread $SLURM_CPUS_PER_TASK \
  ${gene}_pre_align.fasta \
  > ${ALIGNED_DIR}/${gene}_align.fasta \
  2> ${LOG_DIR}/log_${gene}-mafft.out

cd $CLEANED_ALIGNED_DIR

java -jar $APPS/phyutility/phyutility.jar \
  -clean 0.10 \
  -in ${ALIGNED_DIR}/${gene}_align.fasta \
  -out ./tmp/${gene}_cleaned.fasta

sed -i 's/>_R_/>/g' ./tmp/${gene}_cleaned.fasta

python3 $SCRIPT_PATH/sequence_handler.py \
  -f tmp/${gene}_cleaned.fasta \
  --in_type gene \
  --by_gene \
  --min_length ${MIM_SEQ_LENGTH}
