#!/usr/bin/env bash

#SBATCH --job-name=main
#SBATCH --partition=cpuq
#SBATCH --mail-type=NONE
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --chdir=/scratch/alessandro.vinceti/Reduced_Templates

set -o errexit
set -o nounset
set -o pipefail

echo 'Download sgRNA-level fold-change datasets'
snakemake -s workflows/data_download/Snakefile --profile config/config.yml
echo 'Finished downloading data'
