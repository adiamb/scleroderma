args = commandArgs(trailingOnly=TRUE)
file = args[2]
name = args[1]


package_path = "/path/to/packages/"
library("ade4",lib.loc=package_path)
library("seqinr",lib.loc=package_path)

   
data = read.fasta(file)

# Trim subisotype designations:
names(data) = gsub("IgG1_1","IgG1",names(data)) 
names(data) = gsub("IgG1_2","IgG1",names(data)) 
names(data) = gsub("IgG2_1","IgG2",names(data)) 
names(data) = gsub("IgG2_2","IgG2",names(data)) 
names(data) = gsub("IgG3_1","IgG3",names(data)) 
names(data) = gsub("IgG3_2","IgG3",names(data)) 
names(data) = gsub("IgG4_1","IgG4",names(data)) 
names(data) = gsub("IgG4_2","IgG4",names(data)) 
data = lapply(data,function(s) replace(s,s=="-","n"))

# Write results:
n_seqs = length(data)
max_size = 5e5
if (n_seqs > max_size) {
 n_chunks = ceiling(n_seqs/max_size)
 for (i in 1:n_chunks) {
   start = (i-1)*max_size+1
   end = min(i*max_size,n_seqs)
   part_data = data[start:end]
   write.fasta(sequences=part_data, names=names(part_data), nbchar=80,
               file.out=paste(name,"_minCONSCOUNT2_part",i,".fasta",sep=""))
 }
} else {
 write.fasta(sequences=data, names=names(data), nbchar=80,
               file.out=paste(name,"_minCONSCOUNT2.fasta",sep=""))
}