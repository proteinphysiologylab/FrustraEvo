args = commandArgs(trailingOnly=TRUE)

job_id = args[1]
mode= args[2]
PDBStructure=args[3]

file_ic<-paste0("static/jobs/",job_id,"/pdb/FrustraEvo_",job_id,"/OutPutFiles/IC_",mode,"_",job_id)

ICData<-read.table(file_ic, header = T, stringsAsFactors = F)

OutputJML<-paste0("static/jobs/",job_id,"/pdb/FrustraEvo_",job_id,"/OutPutFiles/IC_",mode,"_",job_id,".jml")

blacks<-ICData[which(as.numeric(ICData$ICFrust)<=0.5), "Res"]
blacks<-paste(blacks, sep=",", collapse=",")

#load_pdb<-paste0("background white; load static/jobs/20221215113249687405/pdb/FrustraEvo_",job_id,"/Data/",PDBStructure,".done/VisualizationScrips/",PDBStructure,".pdb;  select all; cartoon only;")
write.table(paste0("background white; load static/jobs/",job_id,"/pdb/FrustraEvo_",job_id,"/Data/",PDBStructure,".done/VisualizationScrips/",PDBStructure,".pdb"), file=OutputJML, quote = F, row.names = F, col.names = F, sep="")

write.table(paste("select", blacks,";"), file=OutputJML, quote = F, row.names = F, col.names = F, append = T)
write.table(paste("color black;"), file=OutputJML, quote = F, row.names = F, col.names = F, append = T)

red<-ICData[which(as.numeric(ICData$ICFrust>0.5) & ICData$FrustrationState=="MAX"), "Res"]

if(length(red)>0)
{
red<-paste(red, sep=",", collapse=",")
write.table(paste("select", red,";"), file=OutputJML, quote = F, row.names = F, col.names = F, append = T)
write.table(paste("color red;"), file=OutputJML, quote = F, row.names = F, col.names = F, append = T)
}


green<-ICData[which(as.numeric(ICData$ICFrust>0.5) & ICData$FrustrationState=="MIN"), "Res"]
if(length(green)>0)
{
green<-paste(green, sep=",", collapse=",")
write.table(paste("select", green,";"), file=OutputJML, quote = F, row.names = F, col.names = F, append = T)
write.table(paste("color green;"), file=OutputJML, quote = F, row.names = F, col.names = F, append = T)
}

gray<-ICData[which(as.numeric(ICData$ICFrust>0.5) & ICData$FrustrationState=="NEU"), "Res"]
if(length(gray)>0)
{
gray<-paste(gray, sep=",", collapse=",")
write.table(paste("select", gray,";"), file=OutputJML, quote = F, row.names = F, col.names = F, append = T)
write.table(paste("color gray;"), file=OutputJML, quote = F, row.names = F, col.names = F, append = T)
}
