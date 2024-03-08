#! /bin/bash

mkdir -p stats/tmp
file_suffix=$1

for file in *$file_suffix
  do echo $file
  gene=$(echo ${file/$file_suffix})
  echo $gene
  AMAS.py summary -i $file -o stats/$gene.txt -f fasta -d dna
done

cd stats
for file in *.txt
  do awk NR==2 $file > tmp/$file
done

ls *.txt | head -n 1 | while read file
 do awk NR==1 $file > heading.txt
done

cat heading.txt tmp/* > final.txt
