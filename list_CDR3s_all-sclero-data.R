args = commandArgs(trailingOnly=TRUE)
files = args

package_path = '/path/to/packages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)

participants = #######

AA_CDR3s <- sapply(participants,function(x) NULL)
for (file in files) {  
  participant = strsplit(file,"_")[[1]][1]
  data = readChangeoDb(file)
  AA_CDR3s[[participant]] = c(AA_CDR3s[[participant]],unique(data$AA_CDR3))
}
all_AA_CDR3s = unique(unlist(AA_CDR3s))

### Write result file:
saveRDS(AA_CDR3s,"AA_CDR3s.RDS")
saveRDS(all_AA_CDR3s,"all_AA_CDR3s.RDS")