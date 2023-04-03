#!/bin/bash
#SBATCH -o ./debug/output_%A.txt
#SBATCH -e ./debug/errors_%A.txt
#SBATCH -J blast
#SBATCH --qos=debug
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6

module load gcc/10.1.0 blast/2.13.0 python

blast_file="$1"
fasta_refseq="$2"
OUT=output/filtered_results

# 1st filtering
python filter_seqs.py $blast_file $OUT 

# get full sequences
db="/apps/BLAST/2.6.0/db_25_03_22/nr"
blastdbcmd -db $db -entry_batch $OUT/${blast_file}.protcoverage70.filtered.seqnames -out $OUT/${blast_file}.protcoverage70.filtered.fa

#2nd filtering
python filter_seqs2.py  $OUT/${blast_file}.protcoverage70.filtered.fa

# add refseq
cat $fasta_refseq $OUT/${blast_file}.protcoverage70.filtered.fa.nosynthetic > $OUT/${blast_file}.protcoverage70.filtered.fa.nosynthetic.withref
