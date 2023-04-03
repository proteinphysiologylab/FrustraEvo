library(data.table)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(stringr)
library(ggh4x)
library(dplyr)
library(seqinr)

#### Read models mapping
mapping_file = fread(
  "../SARS2_FrustraEvo_input_and_output/NSP3_cd21732_betaCoV_PLPro/NSP3_cd21732_betaCoV_PLPro_mapped_positions_unique_MODELS.csv"
)
mapping_file$id = str_remove(mapping_file$PDB_code, '\\.pdb_')


#### Read PDB mapping
mapping_file_pdb = fread(
  "../SARS2_FrustraEvo_input_and_output/NSP3_cd21732_betaCoV_PLPro/NSP3_cd21732_betaCoV_PLPro_mapped_positions_unique_PDB.csv"
)
mapping_file_pdb$id = str_remove(mapping_file_pdb$PDB_code, '\\.pdb_')
mapping_file_pdb

cl1 = Sys.glob(
  "../SARS2_FrustraEvo_input_and_output/NSP3_cd21732_betaCoV_PLPro/frustraevo_output/*/OutPutFilesNSP3_cd21732_betaCoV_PLPro_*/Frustration/*.done/FrustrationData/*.pdb_singleresidue"
)

cl1_f = list()
for (i in cl1) {
  id = str_remove(basename(dirname(dirname(i))), '\\.done')
  cl = str_remove(basename(dirname(dirname(
    dirname(dirname(dirname(i)))
  ))), 'NSP3_cd21732_betaCoV_PLPro_cluster')
  dt = fread(i)
  dt$cll = cl
  if (cl == "allsars2") {
    idd = "nsp3_protein_SARS-CoV-2"
    dt$id = idd
    if (id == "6xa9_A") {
      dt$id2 = idd
    } else{
      dt$id2 = id
    }
    
  } else {
    dt$id = id
    dt$id2 = id
  }
  
  cl1_f [[i]] = dt
}
files = rbindlist(cl1_f)
colnames(files)[1] = "PDB_3D_position"

m2a = merge(subset(files, cll == "allsars2"), mapping_file_pdb)
m2b = merge(subset(files, cll != "allsars2"), mapping_file)
m2 = rbind(m2a, m2b)


m2 = m2 %>% group_by (cll, Subalignment_position) %>%
  mutate(median_FrstIndex = median(FrstIndex),
         mean_FrstIndex = mean(FrstIndex))

mref = subset(
  m2,
  IsRef == "yes" &
    id2 %in% c(
      "AVP25405.1" ,
      "QBF80586.1" ,
      "YP_009555248.1",
      "nsp3_protein_SARS-CoV-2"
    )
)
mref = subset(mref,
              (cll == "all" & id2 == "nsp3_protein_SARS-CoV-2") | cll != "all")


icfrst = Sys.glob(
  "../SARS2_FrustraEvo_input_and_output/NSP3_cd21732_betaCoV_PLPro/frustraevo_output/*/OutPutFilesNSP3_cd21732_betaCoV_PLPro_*/Equivalences/IC_seqfrust.csv"
)
icfrst_f = list()
for (i in icfrst) {
  cl = str_remove(basename(dirname(dirname(i))),
                  'OutPutFilesNSP3_cd21732_betaCoV_PLPro_cluster')
  dt = fread(i, dec = ",")
  
  dt$cll = cl
  
  if (cl == "allsars2") {
    dt$Trimmed_protein_position = dt$Res
  } else{
    dt$Trimmed_protein_position = dt$Res
  }
  icfrst_f [[i]] = dt
}
files2 = rbindlist(icfrst_f, fill = T)

ictot_icseq = merge(mref, files2, by = c("Trimmed_protein_position", "cll"))

##################################
## consensus seq
####

msas = Sys.glob("../msas_PLPro/*.fasta")

res_list = list()
for (n in 1:length(msas)) {
  msa = read.alignment(msas[n],
                       format = "FASTA",
                       forceToLower = F,
                       seqtype = "AA")
  c = consensus(msa)
  cl =  str_remove(str_split_fixed(basename(msas[n]),  '_', 4)[, 4], '.fasta')
  con = data.frame(consensus = c,
                   Subalignment_position = 1:length(c))
  
  if (cl == "all") {
    sub_ictot_icseq = subset(ictot_icseq,  cll == cl)
  } else{
    sub_ictot_icseq = subset(ictot_icseq, Cluster == cl & cll == cl)
  }
  sub_ictot_icseq = merge(sub_ictot_icseq, con, by = "Subalignment_position")
  res_list[[n]] = sub_ictot_icseq
}

nb = subset(ictot_icseq, cll == 'allsars2')
nb$consensus = nb$AA

res_list[[6]] = nb
ictot_icseq = rbindlist(res_list, use.names = T)



ictot_icseq$cll[ictot_icseq$cll == "1"] = "Sarbecovirus"
ictot_icseq$cll[ictot_icseq$cll == "2"] = "Nobecovirus"
ictot_icseq$cll[ictot_icseq$cll == "3"] = "Merbecovirus"
ictot_icseq$cll[ictot_icseq$cll == "4"] = "Embecovirus"
ictot_icseq$cll[ictot_icseq$cll == 'allsars2'] = "SARS2"

ictot_icseq$cll = factor (
  ictot_icseq$cll,
  levels = c(
    "Embecovirus",
    "Merbecovirus",
    "Nobecovirus",
    "Sarbecovirus",
    "SARS2"
  )
)
ictot_icseq$threshold = ifelse(ictot_icseq$ICTot > 0.5, 'ICTot > 0.5', 'ICTot <= 0.5')
r = unique(subset(ictot_icseq, cll == "SARS2")[, c("Subalignment_position", "PDB_3D_position")])
colnames(r)[2] = "SARS2_ref"

ictot_icseq = merge(ictot_icseq, r)
ictot_icseq$median_FrstIndex2  = ictot_icseq$median_FrstIndex
ictot_icseq[ictot_icseq$ICTot <= 0.5, ]$median_FrstIndex2 = NA




##### per region of interest
dom = data.frame(
  domain = c(
    rep("Ub binding site S2", 15),
    rep("Ub binding site S1", 24) ,
    rep("Zinc ion binding site", 4),
    rep("Catalytic site", 5)
  ),
  #,
  # rep("Inhibition loop", 6)),
  SARS2_ref = c(
    62,
    65,
    66,
    69,
    70,
    73,
    75,
    77,
    128,
    175,
    177,
    179,
    180,
    200,
    203,
    157,
    164,
    166,
    167,
    170,
    171,
    174,
    185,
    187,
    197,
    198,
    199,
    203,
    206:209,
    222:225,
    232,
    247,
    269,
    189,
    192,
    224,
    226,
    111,
    272,
    286,
    106,
    109#,
    #266:271
  )
)




ictot_icseq2 = merge(ictot_icseq, dom, by = "SARS2_ref")

ictot_icseq2$SARS2_ref = as.character(ictot_icseq2$SARS2_ref)
ictot_icseq2$SARS2_ref = factor(ictot_icseq2$SARS2_ref, levels = as.character(sort(as.numeric(
  unique(ictot_icseq2$SARS2_ref)
))))


s = subset(ictot_icseq, SARS2_ref %in% dom$SARS2_ref)

ictot_icseq2$ICseq = as.numeric(str_replace(ictot_icseq2$ICseq, ',', '\\.'))

g1 = ggplot(subset(ictot_icseq2, cll != "All clusters"),
            aes(x = SARS2_ref, y = cll)) +
  geom_tile(
    aes(fill = median_FrstIndex2),
    linejoin = "round",
    width = 0.96,
    height = 0.96,
    color = 'grey20'
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_text(aes(label  = consensus, size = ICseq), fontface = "bold") +
  
  coord_cartesian(clip = "off") +
  scale_fill_gradientn(
    values = c(1, .79, .5, 0.11, 0.0909, 0),
    colours = c("green", "grey", "grey", "grey", "red", "red"),
    name = "SRLF median index",
    limits = c(-1.2, 1.1),
    oob = scales::squish,
    na.value = "white"
  ) +
  scale_color_gradient2(
    mid = "grey",
    high = "black",
    midpoint = 0.5,
    limits = c(0, log2(20)),
    oob = scales::squish
  ) +
  theme_bw() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "top") +
  facet_grid( ~ domain, scales = "free_x") +
  force_panelsizes(cols = c(0.09, 0.45, 0.25, 0.09))


saveRDS(g1, 'figure5C.RDS')
