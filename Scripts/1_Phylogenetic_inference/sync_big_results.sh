#!/bin/bash
#BATCH --workdir=/home/azu10kg
#SBATCH --job-name=Sync
#SBATCH --partition=medium
#SBATCH --output=sync_logs.txt
#SBATCH --error=sync_logs.txt
#SBATCH --mem=10M
#SBATCH --open-mode=append
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

now() { echo "$(date +"%F %T"): $1 "; }

echo "###############################"
echo ""
now "Beginning sync"
echo ""

$APPS/rclone/./rclone sync ~/results/  OneDrive:Gruffalo -v

echo ""
now "Sync completed"

echo ""
echo "###############################"
echo ""
