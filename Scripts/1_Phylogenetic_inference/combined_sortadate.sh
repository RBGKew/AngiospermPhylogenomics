#! /bin/bash


source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh
source /home/${USER}/.bashrc
conda activate sortadate

outgroup=$(awk 'NR==1{print;exit}' ${SCRIPT_PATH}/outgroup.txt)
name="combined"
trees="./trees"
flend=".tre"

sh ${SCRIPT_PATH}/combined_aggregate_treeshrink_output.sh $SORTADATE_DIR

cp ${ASTRAL_DIR}/astral.tre ${SORTADATE_DIR}/astral.tre

mkdir -p  $SORTADATE_DIR
cd $SORTADATE_DIR

Rscript ${SCRIPT_PATH}/Root_gene_trees_sortadate.R astral.tre $trees $outgroup
sed "s/'\[[^][]*\]'//g" spt_rooted.tre > clean_tree.tre

echo "Get var length"
python ${SORTADATE_PATH}/src/get_var_length.py $trees/rooted --flend $flend --outf ${name}_var
head -n 5 ${name}_var

echo "Get Bipartitions"
python ${SORTADATE_PATH}/src/get_bp_genetrees.py $trees/rooted clean_tree.tre --flend $flend --outf ${name}_bp.txt
head -n 5 ${name}_bp.txt

python ${SORTADATE_PATH}/src/get_genetree_clades.py $trees/rooted clean_tree.tre --flend $flend --outf ${name}_bp_clades

echo "Combining results"
python ${SORTADATE_PATH}/src/combine_results.py ${name}_var ${name}_bp.txt --outf ${name}_comb.txt
head -n 5 ${name}_comb.txt

n_trees=$(awk '{x++} END {print x}' ${name}_var)

echo $n_trees
echo "Sorting results"
python ${SORTADATE_PATH}/src/get_good_genes.py --order 3,1,2 --max $n_trees --outf ${name}_gg.txt ${name}_comb.txt
head -n 50 ${name}_gg.txt

