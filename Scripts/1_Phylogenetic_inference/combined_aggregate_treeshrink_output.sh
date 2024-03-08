#! /bin/bash

gene_list=$BIGTREE_REPO/gene_list.txt 

treeshrink_dir="07_treeshrink"
genes_dir="genes"
dest_dir=$1

mkdir -p $dest_dir/trees
mkdir -p $dest_dir/aligns

cat $gene_list | while read gene
  do gene_dir=$treeshrink_dir/$genes_dir/$gene
  align_path="$gene_dir/output.fasta"
  tree_path="$gene_dir/output.tre"
  if [ -s "$align_path" ] && [ -s "$tree_path" ]; then
    cp $align_path $dest_dir/aligns/$gene.fasta
    cp $tree_path $dest_dir/trees/$gene.tre
  else
    echo $gene >> dest_dir/missing_genes.txt
  fi
done
