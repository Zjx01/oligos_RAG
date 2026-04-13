#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=cpu
#SBATCH --mem=60G
#SBATCH --output=./message
#SBATCH --time  200:0
#SBATCH --exclusive 
mafft --thread 20  --auto results/rsvb_clean.fasta > results/rsvb_clean.aln.fa 
