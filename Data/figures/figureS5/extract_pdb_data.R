library(data.table)
library(dplyr)
library(stringr)

d = "./sp*.gz"
prot = lapply(Sys.glob(d),fread)

prot = rbindlist(prot)
prot100 = subset(prot, Pident >99)
selpdb = subset(prot100, PDB_code %in% c("1xd2.cif", "3tkl.cif",  "6bcb.cif",'7mge.cif'))
mapping_file_motifs_refs  = fread('./mapping_ras_positions_refs.txt')


selpdb$id = str_split_fixed(selpdb$Protein_accession, "\\|", 3)[,2]
mapping_file_motifs_refs2 = left_join(mapping_file_motifs_refs, selpdb)

fwrite(unique(mapping_file_motifs_refs2[,c("id", "RAS_ref_position", "Alignment_position", "aa",
                                    "Protein_position", "AA_res", "domain", "sdp","PDB_code",
                                    "PDB_chain", "PDB_3D_position", "PDB_aa")]), 'mapping_ras_positions_pdb_motifs.txt')


prot100$id1 = str_split_fixed(prot100$Protein_accession, "\\|", 3)[,2]


contacts = fread("/media/victoria/VICTORIA_EXTERNAL_DISK/CORONAVIRUS/frustration/paper_figs/ras/network_contacts/RAS_contacts_motifs_sdps_ALL.txt")
contacts


subprot = prot100[, c("id1", "Protein_position", "PDB_3D_position", "PDB_code","PDB_chain",
                    "PDB_aa", "Interaction_type", "PDB_interacting_chain", "PDB_interacting_3D_position")]
subprot1 = subprot
subprot2 = subprot

colnames(subprot1) = c("id1", "Res1", "PDB_3D_position1", "PDB_code",
                       "PDB_chain1",  "PDB_aa1", "Interaction_type1",
                       "PDB_interacting_chain1", "PDB_interacting_3D_position1")
colnames(subprot2) = c("id2", "Res2", "PDB_3D_position2", "PDB_code",
                       "PDB_chain2",  "PDB_aa2", "Interaction_type2",
                       "PDB_interacting_chain2", "PDB_interacting_3D_position2")

library(dplyr)
contacts2 = unique(left_join(left_join(contacts, subprot1), subprot2))
contacts2

contacts_selected = subset(contacts2, PDB_code %in% c("1xd2.cif", "3tkl.cif", '7mge.cif',"6bcb.cif"))
fwrite(contacts_selected, './RAS_contacts_motifs_sdps_ALL_annotatedPDB_filtered.txt.gz',
       compress = "gzip")

