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
armA_pats = ###############
armB_pats = #################
healthy_pats = #############
id_to_dummy =  c(paste("SR",1:length(armA_pats),sep=""),
                 paste("SP",1:length(armB_pats),sep=""),
                 paste("H",1:length(healthy_pats),sep=""))
names(id_to_dummy) = c(armA_pats,armB_pats,healthy_pats)

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

##### Plots
### Mutation usage time courses
iso = "all"
mutation_usage_timecourse = readRDS(paste(path,"mutation_usage_timecourse_",iso,".RDS",sep=""))
nmut_max = 10
nmuts=as.character(0:nmut_max)
nmut_colors = rev(brewer.pal(11,"RdYlBu"))
names(nmut_colors) = as.character(nmuts)

for (pat in sort(names(mutation_usage_timecourse))) {
  if (save_eps) {
    pdf(file=paste(fig_path,"/",id_to_dummy[pat],"_mutation_usage_timecourses.eps",sep=""), 
        width=3.2,height=2,paper="a4",
        colormodel="rgb",pointsize=10)
  }
  par(mar=c(5,6,2,2))
  y=mutation_usage_timecourse[[pat]][nmuts]
  max_y = max(sapply(y,function(s) if (s[1]==0) {NA} else {max(s)/s[1]}),na.rm=T)
  xlabel="Time (weeks)"
  ylabel="Number of sequences\nrelative to baseline"
  xlimit=c(-1,110)
  ylimit=c(0,9)
  logaxis=''
  new_plot=T
  for (i in 1:length(nmuts)) {
    nmut = nmuts[i]
    X = as.numeric(names(y[[nmut]]))
    Y = as.vector(y[[nmut]])
    if (Y[1]==0) {Y=Y/(5e-2+Y[1])} else {Y=Y/Y[1]} # normalize to baseline value
    if (new_plot) {
        plot(X[order(X)],Y[order(X)],xlim=xlimit,
             xlab=xlabel,ylab=ylabel,
             ylim=ylimit,log=logaxis,
             axes=F,ann=F,type='o',
             pch=20,col=nmut_colors[nmut])
        new_plot=F
      } else {
        lines(X[order(X)],Y[order(X)],type='o',
              pch=20,col=nmut_colors[nmut])
      }
  }
  if (id_to_dummy[pat] %in% c("SR6")) {xtick_labels=T; xlabel="Time (weeks)"} else {xtick_labels=F; xlabel=""}
  axis(side = 1, at = axTicks(1), labels = xtick_labels)
  axis(side = 2, at = axTicks(2), label = axTicks(2),las=1)
  box(which = "plot", lty = "solid")
  mtext(side = 1, text = xlabel, line = 2)
  abline(1,0,lty=2)
  legend_pos="topright"
  if (id_to_dummy[pat]=="SR6") {legend_pos="topleft"}
  legend(legend_pos,id_to_dummy[pat],bty='n')
  if (save_eps) dev.off()
}

if (save_eps) {
  pdf(file=paste(fig_path,"/mutation_usage_timecourses_legend.eps",sep=""), 
      width=20,height=3,paper="a4",
      colormodel="rgb",pointsize=10)
}
par(mar=c(0,0,0,0))
plot.new()
legend("center",names(nmut_colors),pch=rep(15,length(nmut_colors)),
       col=nmut_colors,bty='n',ncol=6)
if (save_eps) dev.off()