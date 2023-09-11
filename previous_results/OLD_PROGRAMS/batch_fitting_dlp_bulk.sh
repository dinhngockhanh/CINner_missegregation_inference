#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH -A iicd
#SBATCH -J ZIJIN
#       SBATCH -t 5-00:00:00
#SBATCH --cpus-per-task=32
#       SBATCH --mem=700gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zx2406@columbia.edu

module load R

pwd

echo "Launching an R run"
date

R CMD BATCH --no-save --vanilla program_DLPandBULK_Misseg_and_chrs_Zijin.r routput
