#! /bin/bash

# Arguments:
# 1: sorta_date_path
sorta_path=$APPS/SortaDate
# 2: naming
name=$1
# 3: outgroups # If multiple, provided within quotes, comma separated
outgroup=$2
# 4: gene_trees_folder
# 5: species tree
sp_tree=$3

trees=$4

flend=$5

sed "s/'\[[^][]*\]'//g" $sp_tree > clean_tree.tre

Rscript ~/scripts/bigtree/Root_gene_trees_sortadate.R clean_tree.tre $trees $outgroup

echo "Get var length"
python ${sorta_path}/src/get_var_length.py $trees/rooted --flend $flend --outf ${name}_var
head -n 5 ${name}_var

echo "Get Bipartitions"
python ${sorta_path}/src/get_bp_genetrees.py $trees/rooted spt_rooted.tre --flend $flend --outf ${name}_bp.txt
head -n 5 ${name}_bp.txt

python ${sorta_path}/src/get_genetree_clades.py $trees/rooted spt_rooted.tre --flend $flend --outf ${name}_bp_clades

echo "Combining results"
python ${sorta_path}/src/combine_results.py ${name}_var ${name}_bp.txt --outf ${name}_comb.txt
head -n 5 ${name}_comb.txt

n_trees=$(awk '{x++} END {print x}' ${name}_var)

echo $n_trees
echo "Sorting results"
python ${sorta_path}/src/get_good_genes.py --order 3,1,2 --max $n_trees --outf ${name}_gg.txt ${name}_comb.txt
head -n 5 ${name}_gg.txt

