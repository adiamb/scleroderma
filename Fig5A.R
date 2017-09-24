save_eps=T
path=####
fig_path=####
armA_pats = ###############
armB_pats = #################
healthy_pats = #############
id_to_dummy =  c(paste("SR",1:length(armA_pats),sep=""),
                 paste("SP",1:length(armB_pats),sep=""),
                 paste("H",1:length(healthy_pats),sep=""))
names(id_to_dummy) = c(armA_pats,armB_pats,healthy_pats)

### Subisotype
isotype_usage_timecourse = readRDS(paste(path,"subisotype_usage_timecourse.RDS",sep="")) # made using script Fig6A_subisotype-usage-timecourses.R
isotype_cor_timecourse = list()
for (pat in pats) {
  df = data.frame()
  weeks = sort(names(isotype_usage_timecourse[[pat]][[1]]))
  for (week in weeks) {
    isos = sort(names(isotype_usage_timecourse[[pat]]))
    for (iso in isos) {
      df[week,iso] = isotype_usage_timecourse[[pat]][[iso]][week]
    }
  }
  df[is.na(df)] = 0
  for (week in weeks) {
    isotype_cor_timecourse[[pat]][week] = cor(as.numeric(df[weeks[1],]),as.numeric(df[week,]))
  }
}

### VDJ
iso = "all"
gene_usage_timecourse = readRDS(paste(path,"VDJcomb_usage_timecourse_",iso,".RDS",sep="")) # made using script VDJ_comb_usage_timecourse.R
gene_cor_timecourse = list()
for (pat in pats) {
  df = data.frame()
  weeks = sort(names(gene_usage_timecourse[[pat]][[1]]))
  for (week in weeks) {
    genes = sort(names(gene_usage_timecourse[[pat]]))
    for (gene in genes) {
      df[week,gene] = gene_usage_timecourse[[pat]][[gene]][week]
    }
  }
  df[is.na(df)] = 0
  for (week in weeks) {
    gene_cor_timecourse[[pat]][week] = cor(as.numeric(df[weeks[1],]),as.numeric(df[week,]))
  }
}
### Mutation number
iso = "all"
mutation_usage_timecourse = readRDS(paste(path,"mutation_usage_timecourse_",iso,".RDS",sep="")) # made using script Fig6B_mutation-usage-timcourses.R
mutation_cor_timecourse = list()
for (pat in pats) {
  df = data.frame()
  weeks = sort(names(mutation_usage_timecourse[[pat]][[1]]))
  for (week in weeks) {
    nmuts = as.character(0:10)
    for (nmut in nmuts) {
      df[week,nmut] = mutation_usage_timecourse[[pat]][[nmut]][week]
    }
  }
  df[is.na(df)] = 0
  for (week in weeks) {
    mutation_cor_timecourse[[pat]][week] = cor(as.numeric(df[weeks[1],]),as.numeric(df[week,]))
  }
}
### Mutation identity
iso = "all"
mutationIDs_byWeek = readRDS(paste(path,"mutationIDs_byWeek_",iso,".RDS",sep="")) # made using script mutation-identities-timecourse.R
frac_mutIDs_rel_bl = list()
for (pat in pats) {
  weeks = sort(names(mutationIDs_byWeek[[pat]]))
  mutIDs_bl = mutationIDs_byWeek[[pat]][[weeks[1]]]; mutIDs_bl = mutIDs_bl[!grepl("__$",mutIDs_bl)]
  for (week in weeks) {
    mutIDs = mutationIDs_byWeek[[pat]][[week]]; mutIDs = mutIDs[!grepl("__$",mutIDs)]
    frac_mutIDs_rel_bl[[pat]][week] = length(intersect(mutIDs,mutIDs_bl))/length(mutIDs_bl)
  }
}


### Heatmap (Fig. 5A)
isotype_use_corr=unlist(sapply(isotype_cor_timecourse[pats], function(s) {
  names(s)=as.character(round(as.numeric(names(s))))
  s=s[order(as.numeric(names(s)))]
  return(s)}))
VDJ_use_corr=unlist(sapply(gene_cor_timecourse[pats], function(s) {
    names(s)=as.character(round(as.numeric(names(s))))
    s=s[order(as.numeric(names(s)))]
    return(s)}))
mutation_hist_corr=unlist(sapply(mutation_cor_timecourse[pats], function(s) {
    names(s)=as.character(round(as.numeric(names(s))))
    s=s[order(as.numeric(names(s)))]
    return(s)}))
mutation_frac_bl=unlist(sapply(frac_mutIDs_rel_bl[pats], function(s) {
    names(s)=as.character(round(as.numeric(names(s))))
    s=s[order(as.numeric(names(s)))]
    return(s)}))
df = data.frame(isotype_use_corr,VDJ_use_corr,mutation_hist_corr,mutation_frac_bl)

annot = data.frame(row.names=rownames(df),
                   Week=as.numeric(sapply(strsplit(rownames(df),"\\."),function(s) s[2])),
                   Participant=id_to_dummy[sapply(strsplit(rownames(df),"\\."),function(s) s[1])])
patient_colors = brewer.pal(length(unique(annot$Participant)),"Dark2"); names(patient_colors)=sort(unique(annot$Participant))
annotation_colors = list(Participant=patient_colors,Week=c("gray90","gray10"))
variable_descr = c(isotype_use_corr="Isotypes (corr. bl.)",
                   VDJ_use_corr="VDJ combinations (corr. bl.)",
                   mutation_hist_corr="Mutation loads (corr. bl.)",
                   mutation_frac_bl="Specific mutations (frac. bl.)")

if (save_eps) {
  pdf(file=paste(fig_path,"corr-to-bl-heatmap.eps",sep=""), 
      width=20,height=6,
      paper="a4",onefile=F,
      colormodel="rgb",pointsize=10)
}
pheatmap(t(df),annotation_col=annot,cluster_rows=F,cluster_cols=F,
         gaps_col=cumsum(table(annot$Participant)),gaps_row=1:3,
         show_colnames=F,annotation_colors=annotation_colors,
         labels_row=variable_descr[colnames(df)],
         border_color=NA,cellwidth=5,cellheight=8,fontsize=10,fontsize_number=10)
if (save_eps) dev.off()