library('RWebLogo')
N=read.delim("long.txt", stringsAsFactors=F, header=F)
aln=readLines('Logo.fasta')
# Generate the logo in the file
weblogo(seqs=aln, format='png',stacks.per.line = N$V1, color.scheme = 'chemistry', file.out='seqlogo.png')
dev.off()