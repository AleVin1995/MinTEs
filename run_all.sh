#!/usr/bin/env bash

#SBATCH --job-name=main
#SBATCH --partition=cpuq
#SBATCH --mail-type=NONE
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --chdir=/scratch/alessandro.vinceti/Reduced_Templates

source ~/.bashrc
conda activate MinTEs

echo 'Download sgRNA-level fold-change datasets'
snakemake -s workflows/data_download/Snakefile --profile profile
echo 'Finished downloading data'

echo 'Data preprocessing'
for dataset in "Achilles" "Score"
do
    echo $dataset > project.txt
    snakemake -s workflows/data_preprocessing/Snakefile --profile profile

    N_datasets=$(ls resources/BF/*/Project_"$dataset"_BF.tsv | wc -l)

    while [[ $N_datasets -lt 2 ]]
    do
        sleep 30
        N_datasets=$(ls resources/BF/*/Project_"$dataset"_BF.tsv | wc -l)
    done
done

rm project.txt
echo 'Assembling and scaling of sgRNA-level Bayes factor datasets'

rm slurm*
rm -r .snakemake/
rm -r log/