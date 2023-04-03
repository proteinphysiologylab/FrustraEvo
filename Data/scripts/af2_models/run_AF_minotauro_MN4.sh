#!/bin/bash
#SBATCH --job-name=AlphaFold
#SBATCH --output=debug/ColabFold_%j.out
#SBATCH --error=debug/ColabFold_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH -D .
#SBATCH --qos=bsc_ls

module purge && module load singularity alphafold

Path=$(pwd) 

fasta="$1"
output="$2"

echo "Start at $(date)"
echo "-------------------------"


bsc_alphafold --fasta_paths=$fasta --output_dir=$output --max_template_date=2022-03-01


echo "End at $(date)"

