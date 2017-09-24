library(alakazam)
extract_Vgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  s=sapply(w,function(x) {paste(strsplit(x,"-")[[1]][1:2],collapse="-")})
  return(as.vector(s))
}

path = ####
fig_path = ####
save_eps = T
calculate = T

sclero_pats = ####
iso = "all"
if (calculate) {
  mutationIDs_byWeek = list()
  for (pat in sclero_pats) {
    file = paste(path,"/",pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-change-table.tab",sep="")
    data0 = readChangeoDb(file)
    data0$ISOTYPE = gsub("1|2|3|4","",data0$PRIMER)
    if (iso != "all") data0 = subset(data0,ISOTYPE==iso)
    data0$Vgene = extract_Vgene(data0$V_CALL_GENOTYPED)
    mutationIDs_byWeek[[pat]] = list()
    
    for (week in sort(unique(data0$WEEK_CORRECTED))) {
    data = subset(data0,WEEK_CORRECTED==week)
    AA_CHANGES_byVsegment = sapply(unique(data$Vgene),function(s) {
      mutations=unlist(strsplit(data$AA_CHANGES[data$Vgene==s],","))
      mutations=mutations[!is.na(mutations) & mutations!="NA" & mutations!=""]
      mutations=paste(s,mutations,sep="__")
      return(mutations)
    })
    AA_CHANGES_byVsegment = unique(unlist(AA_CHANGES_byVsegment))
    mutationIDs_byWeek[[pat]][[as.character(week)]] = AA_CHANGES_byVsegment
    }
    saveRDS(mutationIDs_byWeek,paste(path,"mutationIDs_byWeek_",iso,".RDS",sep=""))
  }
}