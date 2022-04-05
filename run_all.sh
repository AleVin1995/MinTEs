#!/usr/bin/env bash

#SBATCH --job-name=main
#SBATCH --partition=cpuq
#SBATCH --mail-type=NONE
#SBATCH --time=02:00:00
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
        N_jobs_assemble=0

        while [[ $N_datasets -lt 2 && $N_jobs_assemble -eq 0 ]]
        do
            sleep 30
            N_datasets=$(ls resources/BF/*/ | grep $dataset | grep tsv | wc -l)
            N_jobs_assemble=$(squeue | grep assemble | wc -l)
        done

        sleep 30
        N_cells_BF=$(head -1 resources/BF/*/Project_"$dataset"_BF.tsv | wc -w)
        
        if [[ $N_cells_FC -gt $N_cells_BF ]]
        then
            rm resources/BF/*/Project_"$dataset"_BF.tsv
        fi
    done
done

rm project.txt
echo 'Assembling and scaling of sgRNA-level Bayes factor datasets'

#rm slurm*
#rm -r .snakemake/
#rm -r log/