#!/bin/bash
#SBATCH --job-name=FrstEvo
#SBATCH --output=/gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/debug/FrustraEvo_%j.out
#SBATCH --error=/gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/debug/FrustraEvo_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH -D .
#SBATCH --qos=debug

# load necessary modules
module load mkl python

# activate virtual environment
source /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/venv_frustraevo/bin/activate

module load R modeller hmmer

# set python library path
PYTHONPATH=/gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/venv_frustraevo/lib

# input
jobID="$1"
path_to_r=/gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/FrustraEvo2.0/Scripts
msa_file="$2"
list_ids="$3"
seqrefid="$4"
pdbs_path="$5"

# run frustraEvo
echo "Start at $(date)"
echo "-------------------------"

python /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/FrustraEvo2.0/run_logo.py $jobID $path_to_r $msa_file $list_ids $seqrefid $pdbs_path

echo "End at $(date)"

