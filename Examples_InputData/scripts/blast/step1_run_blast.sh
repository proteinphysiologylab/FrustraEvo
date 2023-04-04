#!/bin/bash
#SBATCH -o ./debug/output_%A.txt
#SBATCH -e ./debug/errors_%A.txt
#SBATCH -J blast
#SBATCH --qos=debug
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6

module load gcc/10.1.0 blast/2.13.0

export BLASTDB_LMDB_MAP_SIZE=100000000

db="/apps/BLAST/2.6.0/db_25_03_22/nr"


blastp -query "$1" -out "$2" -max_target_seqs 10000 -evalue 0.05 -num_threads 4 -db $db -outfmt "6 qseqid sseqid staxid evalue pident bitscore qcovs sallseqid qlen slen qstart qend sstart send length nident qseq sseq" 
