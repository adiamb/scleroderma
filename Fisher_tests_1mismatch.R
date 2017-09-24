args = commandArgs(trailingOnly=TRUE)
START = as.numeric(args[1])
BIN_SIZE = as.numeric(args[2])
END = min(START+BIN_SIZE-1,2410035)
outfile = args[3]
AA_CDR3s_file = args[4]
AA_CDR3s = readRDS(AA_CDR3s_file)
all_AA_CDR3s_file = args[5]
all_AA_CDR3s = readRDS(all_AA_CDR3s_file)

package_path = '/path/to/packages/'
library('ggplot2'); library('plyr'); library('reshape2'); library('scales')
library('lazyeval',lib.loc=package_path)
library('igraph',lib.loc=package_path)
library('alakazam',lib.loc=package_path)
library('stringdist', lib.loc=package_path)

healthy_participants = #####
sclero_participants = #####
participants = c(healthy_participants,sclero_participants)

SCLERO = as.numeric(participants %in% sclero_participants)
pvals = data.frame()
for (i in START:END) { 
  CDR3 = all_AA_CDR3s[i]
  PRESENT = sapply(participants,function(participant) as.numeric(ain(CDR3,AA_CDR3s[[participant]],method="hamming",maxDist=1)))
  ## Compute Fisher exact test
  if (length(unique(PRESENT))==1) {
    pvals[CDR3,"raw"] = 1
  } else {	
    pvals[CDR3,"raw"] = fisher.test(SCLERO,PRESENT)$p.value
  }
  ## Record participants in which sequence was observed
  for (participant in participants) {
    pvals[CDR3,participant] = PRESENT[participant]
  }
  ## Record incidences in the two groups
  pvals[CDR3,"incidence_sclero"] = sum(PRESENT[sclero_participants])
  pvals[CDR3,"incidence_healthy"] = sum(PRESENT[healthy_participants])
}
saveRDS(pvals,outfile)