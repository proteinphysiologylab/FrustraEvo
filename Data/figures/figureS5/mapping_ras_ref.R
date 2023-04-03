library(data.table)
library(seqinr)
library(stringr)
library(dplyr)

##################################
msa = read.fasta('../figure4/figure4B/ras_human.fasta',
                 forceDNAtolower = F, seqtype = "AA")

names(msa) = str_split_fixed(names(msa),pattern = '\\|',n = 3)[,2]
names(msa) = str_remove(names(msa), '-1')

seqrefs= c("P01112", "P84077" ,"P01112", "P61586", "P62820")

list_df = list()

for (i in 1:length(msa)){
  s = msa[[i]]
  index_a = 0
  index_s = 0
  aa = c()
  seqpos = c()
  alipos = c()

  for(j in 1:length(s)){
    if (s[j] !='-'){
      index_a = index_a + 1
      index_s = index_s + 1
      aa = c(aa, s[j])
      seqpos = c(seqpos, index_s)
      alipos = c(alipos, index_a)
    } else {
      index_a = index_a + 1
    }
  }
  if (names(msa)[i] %in% seqrefs){
    IsRef = "yes"
  } else {
    IsRef = "no"
  }
  
  df = data.frame(aa, Protein_position = seqpos, Alignment_position = alipos, id= names(msa)[i], IsRef )
  list_df[[i]] = df
  
  
}

mapping_file = rbindlist(list_df)


refpos = subset(mapping_file, id == "P01112")
colnames(refpos)[2] = "RAS_ref_position"
refpos$AA_res = paste(refpos$aa, refpos$RAS_ref_position, sep="")
mapping_file = merge(mapping_file, refpos[,c(2:3,6)], by="Alignment_position", allow.cartesian = T)

refs = data.frame(id = seqrefs, id2 = c( "ALL","ARF","RAS","RHO",  "RAB" ) )


i = c(10, 33, 53,65,68,75,81,83,115,145)
j = c(20, 37, 61,65,68,75,81, 83,121, 148)

v = c()
for (k in 1:length(i)){
  v = c(v,i[k]:j[k])
}
v

df = data.frame(res = v,
                domain = c(rep("G1", 11 ),
                           rep("G2",5),
                           rep("G3", 9),
                           rep('-',5),
                           rep("G4",7),
                           rep("G5",4)),
                sdp = c("no","yes",
                        rep ("no",2),"yes",
                        rep ("no",5),"yes",
                        "no","yes",
                        rep ("no",2),"yes",
                        rep ("no",3),"yes",
                        "no","yes","yes",
                        rep ("no",2),rep ("yes",5),
                        rep ("no",6),"yes",
                        rep ("no",4)))

mapping_file_motifs = merge(mapping_file, df, by.x ="RAS_ref_position", by.y= "res")
mapping_file_motifs_refs = subset(mapping_file_motifs, id %in% seqrefs)
fwrite(mapping_file_motifs_refs, 'mapping_ras_positions_refs.txt')