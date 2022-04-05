#!/usr/bin/env bash

#SBATCH --job-name=main
#SBATCH --partition=cpuq
#SBATCH --mail-type=NONE
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --chdir=/scratch/alessandro.vinceti/Reduced_Templates

source ~/.bashrc
conda activate MinTEs

echo 'Download sgRNA-level fold-change datasets'
snakemake -s workflows/data_download/snakefile --profile workflows/data_download/config
echo 'Finished downloading data'

echo 'Data preprocessing'
for dataset in "Achilles" "Score"
do
    N_cells_FC=$((2*$(head -1 resources/FC/Project_"$dataset"_corrected_FC.tsv | wc -w)))
    N_cells_BF=0

    while [[ $N_cells_FC -gt $N_cells_BF ]]
    do
        echo $dataset > project.txt
        snakemake -s workflows/data_preprocessing/snakefile --profile workflows/data_preprocessing/config

        N_datasets=0
        N_jobs_scaling=0

        while [[ $N_datasets -lt 2 || $N_jobs_scaling -ne 0 ]]
        do
            sleep 30
            N_datasets=$(ls resources/BF/*/ | grep $dataset | grep scaled | wc -l)
            N_jobs_scaling=$(squeue | grep scaling | wc -l)
        done

        sleep 30
        N_cells_BF=$(head -1 resources/BF/*/Project_"$dataset"_scaled_BF.tsv | wc -w)
        
        if [[ $N_cells_FC -gt $N_cells_BF ]]
        then
            rm -r resources/FC/*cells
            rm -r resources/BF
        fi
    done
done

rm project.txt
rm -r resources/FC/*cells
rm -r resources/BF/*/*cells
echo 'Finished data preprocessing'

echo 'Reduced templates (RT) optimisation'
snakemake -s workflows/RT_optimisation/snakefile --profile workflows/RT_optimisation/config
echo 'Finished RT optimisation'

rm slurm*
rm -r .snakemake/
rm -r log/