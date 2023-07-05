#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH -A iicd
#SBATCH -J KNDinh
#     SBATCH -t 12:00:00
#SBATCH -t 5-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem = 700gb
#     SBATCH --exclusive
#     SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=knd2127@columbia.edu

module load R

pwd

echo "Launching an R run"
date

R CMD BATCH --no-save --vanilla program_dlp3.r routput
