#!/bin/bash
#
#SBATCH --job-name=Comb_Astral
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=00_logs/Astral.%A.log

echo The script ordinal_astral.sh was called as: 
echo $0 $@

SCRIPT_PATH="~/scripts/bigtree"
START_DIR=$(pwd -LP)
dataset=$(pwd | awk -F/ '{print $(NF)}')
prefix=$(echo $dataset | awk '{print substr($0,0,6)}')

astral_dir="11_astral_gpu"

mkdir -p $astral_dir
cd $astral_dir

ulimit -c unlimited
ulimit -s unlimited
ulimit -H unlimited

astral_mem=$1
astral_path="$APPS/Astral_mp"

java -Xmx$astral_mem -D"java.library.path=$astral_path/lib/" -jar $astral_path/astral.5.15.5.jar \
-t 2 \
-i ../11_astral/${prefix}-BS30.trees \
-o astral.tre \
2> astral.log \
