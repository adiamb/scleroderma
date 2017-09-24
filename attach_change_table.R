args = commandArgs(trailingOnly=TRUE)
wdir = args[1]
file = args[2]

package_path = '/path/to/packages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)
removeBrackets = function(entries) as.numeric(sapply(entries,function(s) strsplit(s," ")[[1]][1]))

data = readChangeoDb(file)
								
### Record mutations in the V segment reported by IMGT:
# Read IMGT output:
patient = strsplit(file,"_")[[1]][1]
IMGTdata = data.frame()
	IMGTfolders = list.files(wdir,patient)
	IMGTfolders = IMGTfolders[grepl("_minCONSCOUNT2",IMGTfolders) & 
								!grepl("txz",IMGTfolders) & 
								!grepl("tab",IMGTfolders) & 
								!grepl("fasta",IMGTfolders) &
								!grepl("extracted",IMGTfolders)]
	for (IMGTfolder in IMGTfolders) {
		week = strsplit(IMGTfolder,"_")[[1]][2]
		IMGTfile = paste(wdir,"/",IMGTfolder,"/7_V-REGION-mutation-and-AA-change-table.txt",sep="")
		newIMGTdata = read.delim(IMGTfile,stringsAsFactors=F)
		newIMGTdata$SEQUENCE_ID = sapply(newIMGTdata$Sequence.ID, function(s) strsplit(s,"_")[[1]][1])
		newIMGTdata$SEQUENCE_VISIT_ID = paste(newIMGTdata$SEQUENCE_ID,week,sep=";")
		IMGTdata = rbind(IMGTdata,newIMGTdata)
	}

# Parse mutation field:
IMGTdata$NT_CHANGES = as.character(sapply(IMGTdata$V.REGION, function(s) {
  mutations = strsplit(s,"\\|")[[1]]
  nt_mutations = sapply(mutations, function(s) strsplit(s,",")[[1]][1])
  nt_mutations = nt_mutations[!grepl('n',nt_mutations)]
  return(paste0(nt_mutations,collapse=","))
}))
IMGTdata$AA_CHANGES = as.character(sapply(IMGTdata$V.REGION, function(s) {
  mutations = strsplit(s,"\\|")[[1]]
  AA_mutations = sapply(mutations, function(s) strsplit(s,",")[[1]][2])
  AA_mutations = AA_mutations[!grepl('X',AA_mutations)]
  return(paste0(AA_mutations,collapse=","))
}))

IMGTdata = IMGTdata[,c("SEQUENCE_VISIT_ID","NT_CHANGES","AA_CHANGES")]

# Bind mutation counts to the changeo database:
data = merge(data,IMGTdata,by="SEQUENCE_VISIT_ID")

### Write result file:
outfile = gsub("_attach-AA-mutationtypes-pass-corr.tab","_attach-change-table.tab",file)
writeChangeoDb(data,outfile)