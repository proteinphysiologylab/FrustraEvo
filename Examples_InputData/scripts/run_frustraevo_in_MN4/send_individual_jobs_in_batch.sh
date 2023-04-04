#!/bin/bash



#################################
# ALL MODELS full alignment
#################################

# run in sbatch each cluster of each protein
for d in $(echo /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/*)
do
# create outdir
out=${d}/frustraevo_output_allmodels
mkdir -p $out	
# initialize cluster count
cl=0
# iterate over subalignment files (i.e.: number of clusters)
for m in $(echo $d/msas_filtered/trimmedtopdb_*.fasta)
do

# create output dir for cluster
mkdir -p ${out}/all
# change directory to generate frustraevo files in the cluster folder
cd ${out}/all
#for subalignment fasta file, extract list of ids								
grep '>' $m | sed  's/>//' > ${out}/all/list_seqs.txt
# define inputs
jobID=$(basename $d)_all		
msa_file=$m
list_ids=${out}/all/list_seqs.txt		
seqrefid=$(awk -F '\t' '{ if($4 =="yes" && $3 == '1') print $2}' ${d}/pdbs_and_models_list.txt)
pdbs_path=${d}/RenamedModels/
# run frustraevo
sbatch /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/run_frustraevo.sh $jobID $msa_file $list_ids $seqrefid $pdbs_path
done
done

#################################
# GRB2 AND PSD95
#################################

# run in sbatch each cluster of each protein
for d in $(echo /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/lehner_paper)
do

for prot in GRB2 PSD95
do
# create outdir
out=${d}/frustraevo_output/${prot}
mkdir -p $out	
# initialize cluster count
#cl=0
# iterate over subalignment files (i.e.: number of clusters)
for m in $(echo $d/msa/filtered_renamed_${prot}.fasta)
do
# update cluster count
#cl=$((cl+1))
# create output dir for cluster
mkdir -p ${out} #/cluster$cl
# change directory to generate frustraevo files in the cluster folder
cd ${out}  #/cluster$cl
#for subalignment fasta file, extract list of ids								
grep '>' $m | sed  's/>//'  | awk -F ' ' '{print $1}' > ${out}/list_seqs.txt

# define inputs
jobID=$prot		
msa_file=$m		
seqrefid=$(head -n 1 ${m} | sed 's/>//' | awk -F ' ' '{print $1}')
pdbs_path=${d}/RenamedModels_${prot}
#basename -s .pdb $(ls ${pdbs_path}) > ${out}/list_seqs.txt
list_ids=${out}/list_seqs.txt

echo $jobID
echo $msa_file
echo $list_ids
echo $seqrefid
echo $pdbs_path
# run frustraevo
sbatch /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/run_frustraevo.sh $jobID $msa_file $list_ids $seqrefid $pdbs_path
done
done
done


#################################
# LEHNER PAPER KRAS NEW PAPER
#################################

# run in sbatch each cluster of each protein
for d in $(echo /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/lehner_paper)
do

for prot in RAS
do
# create outdir
out=${d}/frustraevo_output/${prot}
mkdir -p $out	
# initialize cluster count
#cl=0
# iterate over subalignment files (i.e.: number of clusters)
for m in $(echo /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/ras_family/RAS/msa/filtered_renamed_ras.fasta)
do
# update cluster count
#cl=$((cl+1))
# create output dir for cluster
mkdir -p ${out} #/cluster$cl
# change directory to generate frustraevo files in the cluster folder
cd ${out}  #/cluster$cl
#for subalignment fasta file, extract list of ids								
grep '>' $m | sed  's/>//'  | awk -F ' ' '{print $1}' > ${out}/list_seqs.txt

# define inputs
jobID=$prot		
msa_file=$m		
seqrefid=P01116-2
pdbs_path=/gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/ras_family/RenamedModels
#basename -s .pdb $(ls ${pdbs_path}) > ${out}/list_seqs.txt
list_ids=${out}/list_seqs.txt

echo $jobID
echo $msa_file
echo $list_ids
echo $seqrefid
echo $pdbs_path
# run frustraevo
sbatch /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/run_frustraevo.sh $jobID $msa_file $list_ids $seqrefid $pdbs_path
done
done
done

#################################
# LEHNER PAPER KRAS NEW PAPER ALL MODELS
#################################

# run in sbatch each cluster of each protein
for d in $(echo /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/lehner_paper)
do

for prot in KRAS
do
# create outdir
out=${d}/frustraevo_output/${prot}
mkdir -p $out	
# initialize cluster count
cl=0
# iterate over subalignment files (i.e.: number of clusters)
for m in $(echo $d/msa/filtered_*kras*_*_CAP.fasta)
do
cl=$((cl+1))
# create output dir for cluster
mkdir -p ${out}/cluster$cl
# change directory to generate frustraevo files in the cluster folder
cd ${out}/cluster$cl
#for subalignment fasta file, extract list of ids								
grep '>' $m | sed  's/>//'  | awk -F ' ' '{print $1}' > ${out}/cluster${cl}/list_seqs.txt

# define inputs
jobID=$prot_${cl}		
msa_file=$m		
seqrefid=$(head -n 1 ${m} | sed 's/>//' | awk -F ' ' '{print $1}')
pdbs_path=/gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/lehner_paper/RenamedModels_KRAS
#basename -s .pdb $(ls ${pdbs_path}) > ${out}/list_seqs.txt
list_ids=${out}/cluster${cl}/list_seqs.txt
echo $jobID
echo $msa_file
echo $list_ids
echo $seqrefid
echo $pdbs_path
# run frustraevo
sbatch /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/run_frustraevo.sh $jobID $msa_file $list_ids $seqrefid $pdbs_path

# update cluster count

done

done
done


#################################
# LEHNER PAPER KRAS NEW PAPER ALL MODELS without clusters
#################################

# run in sbatch each cluster of each protein
for d in $(echo /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/lehner_paper)
do

for prot in KRAS
do
# create outdir
out=${d}/frustraevo_output/${prot}
mkdir -p $out	
# initialize cluster count
#cl=0
# iterate over subalignment files (i.e.: number of clusters)
for m in $(echo $d/msa/filtered_*kras*hmm*.fasta)
do
#cl=$((cl+1))
# create output dir for cluster
mkdir -p ${out}#/cluster$cl
# change directory to generate frustraevo files in the cluster folder
cd ${out}#/cluster$cl
#for subalignment fasta file, extract list of ids								
grep '>' $m | sed  's/>//'  | awk -F ' ' '{print $1}' > ${out}/list_seqs.txt #cluster${cl}/list_seqs.txt

# define inputs
jobID=$prot #_${cl}		
msa_file=$m		
seqrefid=$(head -n 1 ${m} | sed 's/>//' | awk -F ' ' '{print $1}')
pdbs_path=/gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/lehner_paper/RenamedModels_KRAS
#basename -s .pdb $(ls ${pdbs_path}) > ${out}/list_seqs.txt
list_ids=${out}/list_seqs.txt #cluster${cl}/list_seqs.txt
echo $jobID
echo $msa_file
echo $list_ids
echo $seqrefid
echo $pdbs_path
# run frustraevo
sbatch /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/run_frustraevo.sh $jobID $msa_file $list_ids $seqrefid $pdbs_path

# update cluster count

done

done
done



#################################
# RAS SUPERFAMILY HUMAN
#################################

prots=(RAB ARF RAS RHO) 
refs=(P62820 P84077 P01112 P61586) 

#for i in ${!prots[*]}; do 
#  echo "${prots[$index]} is in ${refs[$index]}"
#done

for index in ${!prots[*]}
do
d=/gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/ras_family

out=${d}/frustraevo_output/${prots[$index]}
mkdir -p $out

m=$(echo $d/${prots[$index]}/msa/filtered_renamed_*.fasta)
grep '>' $m | sed  's/>//'  | awk -F ' ' '{print $1}' > ${out}/list_seqs.txt
cd ${out}

# define inputs
jobID=${prots[$index]}		
msa_file=$m		
pdbs_path=${d}/RenamedModels
list_ids=${out}/list_seqs.txt
seqrefid=${refs[$index]}

echo $jobID
echo $msa_file
echo $list_ids
echo $seqrefid
echo $pdbs_path

sbatch /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/run_frustraevo.sh $jobID $msa_file $list_ids $seqrefid $pdbs_path

done




#################################
# HEMOGLOBINS alpha beta both
#################################

# run in sbatch each cluster of each protein
for d in $(echo /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/hemoglobins)
do
# create outdir
out=${d}/frustraevo_output
mkdir -p $out	
# initialize cluster count
for cl in Alpha Beta Both
do
# iterate over subalignment files (i.e.: number of clusters)
for m in $(echo $d/msas/${cl}*.fasta)
do

# create output dir for cluster
mkdir -p ${out}/$cl
# change directory to generate frustraevo files in the cluster folder
cd ${out}/$cl
#for subalignment fasta file, extract list of ids								
grep '>' $m | sed  's/>//' > ${out}/${cl}/list_seqs_${cl}.txt
# define inputs
jobID=$(basename $d)_$cl		
msa_file=$m
list_ids=${out}/${cl}/list_seqs_${cl}.txt	

if [ $cl = Alpha ]; then
seqrefid=2dn1-A
else 
seqrefid=2dn1-B
fi
	
pdbs_path=${d}/pdbs/
# run frustraevo
sbatch /gpfs/scratch/bsc08/bsc08453/postdoc/frustration_analysis/frustraevo/run_scripts_sbatch/run_frustraevo.sh $jobID $msa_file $list_ids $seqrefid $pdbs_path
done
done
done



