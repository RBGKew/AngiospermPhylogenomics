#!/bin/bash
#
#SBATCH --job-name=Comb_Astral
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=himem
#SBATCH --output=00_logs/Astral.%A.log
#SBATCH --mem=599G

source ./parameters.sh
source ${SCRIPT_PATH}/combined_standards.sh
source /home/${USER}/.bashrc
conda activate astral_mp

mkdir -p ${ASTRAL_DIR}
cd ${ASTRAL_DIR}

ulimit -c unlimited
ulimit -s unlimited
ulimit -H unlimited

cat ${TREESHRINK_DIR}/genes/*/output.tre > input.trees

Rscript ${SCRIPT_PATH}/root_by_vector.R \
  input.trees \
  ${SCRIPT_PATH}/outgroups.txt

java -Xmx600G \
  -D"java.library.path=$ASTRAL_PATH/lib/" \
  -jar $ASTRAL_PATH/astral.5.15.5.jar \
  -t 2 \
  -C \
  -i rooted_input.trees \
  -o astral.tre \
  2> astral.log

grep "Taxon occupancy: " astral.log > taxon_occupancy.txt

cp astral.tre ~/bigtree_results/bigtree_astral_$(date +"%F").tre

sh ~/scripts/bigtree/sync_big_results.sh
