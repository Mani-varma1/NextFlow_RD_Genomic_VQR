#!/bin/bash
#SBATCH --job-name=NF_UNI
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --output=MainLog.out
#SBATCH --error=MainError.err



nextflow run main.nf -profile slurm


