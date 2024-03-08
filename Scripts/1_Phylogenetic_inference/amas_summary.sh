#!/bin/bash
#SBATCH --job-name=AMAS_Sum
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --mem=500M

### 1: File extension
### 2: Name for summary file

mkdir -p amas_tmp

cp *${1} amas_tmp/

cd amas_tmp

find -empty -delete
AMAS.py summary -f fasta -d dna -i *

cp summary.txt ~/bigtree_results/${2}_amas.txt

sbatch ~/scripts/bigtree/sync_big_results.sh

cd ..

rm -r amas_tmp/

