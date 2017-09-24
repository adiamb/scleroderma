args = commandArgs(trailingOnly=TRUE)
infiles = args[1:(length(args)-1)]
outfile = args[length(args)]

pvals_combined = data.frame()
for (infile in infiles) {
  pvals = readRDS(infile)
  pvals_combined = rbind(pvals_combined,pvals)
}
saveRDS(pvals_combined,outfile)