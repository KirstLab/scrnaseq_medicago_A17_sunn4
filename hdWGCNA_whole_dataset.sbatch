#!/bin/sh
#SBATCH --job-name=hdWGCNA_whole_dataset
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wendellpereira@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300gb
#SBATCH --time=3:50:00
#SBATCH --output=hdWGCNA_whole_dataset_person_%j.out
#SBATCH --qos=kirst

module load R/4.2

Rscript hdWGCNA_whole_dataset.R