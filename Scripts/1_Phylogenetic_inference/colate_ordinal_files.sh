#! /bin/bash

verbose=$1

source ./parameters.sh


aligns_i1="$START_DIR/02_aligned_i1"
aligns_i2="$START_DIR/07_aligned_i2"
trees_i1="$START_DIR/04_genetree_i1"
trees_i2="$START_DIR/09_genetree_i2"
unaligned_seqs="$START_DIR/01_seqs_i1"

colation_log=${PREFIX}_data_colation.log

mkdir -p $FINAL_FILES_DIR/tmp

echo "Data collection started: "$(date '+%d/%m/%Y %H:%M:%S') > $colation_log

cat ~/scripts/bigtree/gene_list.txt | while read gene
  do if [ -f $trees_i2/${gene}_i2.treefile ]; then
    cp $verbose $trees_i2/${gene}_i2.treefile $FINAL_FILES_DIR/${gene}_${PREFIX}_i2.treefile
    echo "tree:"$gene":i2" >> $colation_log
  elif [ -f $trees_i1/${gene}_i1.treefile ]; then
    cp $verbose $trees_i1/${gene}_i1.treefile $FINAL_FILES_DIR/${gene}_${PREFIX}_i1.treefile
    echo "tree:"$gene":i1" >> $colation_log
  else
    echo "tree:"$gene":absent" >> $colation_log
  fi
  if [ -f $aligns_i2/${gene}_i2_aligned.fasta ]; then
    cp $verbose $aligns_i2/${gene}_i2_aligned.fasta $FINAL_FILES_DIR/${gene}_${PREFIX}_i2.fasta
    echo "align:"$gene":i2" >> $colation_log
  elif [ -f $aligns_i1/${gene}_i1_aligned.fasta ]; then
    cp $verbose $aligns_i1/${gene}_i1_aligned.fasta $FINAL_FILES_DIR/${gene}_${PREFIX}_i1.fasta
    echo "align:"$gene":i1" >> $colation_log
  elif [ -f $unaligned_seqs/${gene}.fasta ]; then 
    cp $verbose $unaligned_seqs/${gene}.fasta $FINAL_FILES_DIR/tmp/${gene}.fasta
    echo "align:"$gene":unaligned" >> $colation_log
  else
    echo "align:"$gene":absent" >> $colation_log
  fi
done

if [ "$(ls $DEST_DIR/tmp | wc -l)" -gt 0 ]; then
  echo sequence_handler.py
  python3 $SCRIPT_PATH/sequence_handler.py \
  -f $FINAL_FILES_DIR/tmp/ \
  -v \
  --single_file \
  --in_type gene \
  --out $FINAL_FILES_DIR/${PREFIX}_seqs_not_aligned
fi

rm -r $FINAL_FILES_DIR/tmp/

echo "Data collection ended: "$(date '+%d/%m/%Y %H:%M:%S') >> $colation_log

