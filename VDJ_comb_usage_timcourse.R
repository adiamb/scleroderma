library(RColorBrewer)
library(pheatmap)
library(alakazam)

extract_Vgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  s=sapply(w,function(x) {paste(strsplit(x,"-")[[1]][1:2],collapse="-")})
  return(as.vector(s))
}
extract_Dgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(strsplit(x,"\\*")[[1]][1],"Homsap ")[[1]][2]})
  return(as.vector(w))
}
extract_Jgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(strsplit(x,"\\*")[[1]][1],"Homsap ")[[1]][2]})
  return(as.vector(w))
}



path = ####
fig_path = ####
save_eps = T
calculate = T

sclero_pats = ####
armA_pats = ####
iso = "all"
if (calculate) {
  gene_usage_timecourse = list()
  for (pat in sclero_pats) {
    file = paste(path,"/",pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab",sep="")
    data0 = readChangeoDb(file)
    data0$V_SEGMENT = extract_Vgene(data0$V_CALL_GENOTYPED)
    data0$D_SEGMENT = extract_Dgene(data0$D_CALL)
    data0$J_SEGMENT = extract_Jgene(data0$J_CALL)
    data0 = data0[!is.na(data0$D_SEGMENT),]
    data0$VDJ_COMB = paste(data0$V_SEGMENT,data0$D_SEGMENT,data0$J_SEGMENT,sep="_")
    data0$ISOTYPE = gsub("1|2|3|4","",data0$PRIMER)
    if (iso != "all") data0 = subset(data0,ISOTYPE==iso)
    gene_usage_timecourse[[pat]] = list()
    for (gene in sort(unique(data0$VDJ_COMB))) {
      data1 = subset(data0,VDJ_COMB==gene)
      for (week in sort(unique(data0$WEEK_CORRECTED))) {
        data = subset(data1,WEEK_CORRECTED==week)
        gene_usage_timecourse[[pat]][[gene]][as.character(week)] = nrow(data)
      }
    }
    saveRDS(gene_usage_timecourse,paste(path,"VDJcomb_usage_timecourse_",iso,".RDS",sep=""))
  }
}