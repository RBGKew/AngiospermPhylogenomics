#! /bin/bash
# This script uses two parameters: (1) directory for files (2) file to save
directory=$1
out_file=$2

touch $out_file
for path in $directory/g*_direction.txt
  do file=$(basename $path)
  echo $path
  gene=$(echo ${file/_to_reverse.txt})
  cat $path | while read sample
    do name=$(echo ${sample/">_R_"})
    echo ${name}-$gene >> $out_file
  done
done
