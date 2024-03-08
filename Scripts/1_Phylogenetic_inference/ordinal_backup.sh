#! /bin/bash

START_DIR=$(pwd -LP)
DATASET=$(pwd | awk -F/ '{print $(NF)}')

destination="~/users_area/Bigtree/ordinal/$DATASET/"

mkdir -p $destination

cp *.txt $destination/.

tar -czf $destination/01_seqs_i1.tar.gz 01_seqs_i1/*.fasta

if [ -d 02_aligned_i1 ]; then
 tar -czf $destination/02_aligned_i1.tar.gz 02_aligned_i1/*.fasta
fi

if [ -d 04_genetree_i1 ]; then
 tar -czf $destination/04_genetree_i1.tar.gz 04_genetree_i1/*.treefile
fi

if [ -d 07_aligned_i2 ]; then
 tar -czf $destination/07_aligned_i2.tar.gz 07_aligned_i2/*.fasta
fi

if [ -d 09_genetree_i2 ]; then
 tar -czf $destination/09_genetree_i2.tar.gz 09_genetree_i2/*.treefile
fi

if [ -d 10_trees_aligns ]; then
 tar -czf $destination/10_trees_aligns.tar.gz 10_trees_aligns/*
fi

if [ -d 11_astral ]; then
 cp 11_astral/*.tre $destination/.
fi
