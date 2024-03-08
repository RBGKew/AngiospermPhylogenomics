#! /bin/bash

gene_list=$BIGTREE_REPO/gene_list.txt 
outgroup=$(awk 'NR==1{print;exit}' $BIGTREE_REPO/outgroup.txt)

dest_dir=09_to_calibrate

sh $BIGTREE_REPO/combined_aggregate_treeshrink_output.sh $dest_dir

cp 08_astral/astral.tre $dest_dir/astral.tre

Rscript $BIGTREE_REPO/Root_gene_trees_sortadate.R $dest_dir/astral.tre $dest_dir/trees/ $outgroup
