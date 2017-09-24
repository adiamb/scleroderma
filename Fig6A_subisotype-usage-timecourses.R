library(RColorBrewer)
library(pheatmap)
library(alakazam)
library(tidyr)

sclero_pats = ################
armA_pats = ###############
armB_pats = #################
healthy_pats = #############
id_to_dummy =  c(paste("SR",1:length(armA_pats),sep=""),
                 paste("SP",1:length(armB_pats),sep=""),
                 paste("H",1:length(healthy_pats),sep=""))
names(id_to_dummy) = c(armA_pats,armB_pats,healthy_pats)


path = ####
fig_path = ####
save_eps = T
calculate = T

isotypes = c("IgM","IgA1","IgA2","IgG1","IgG2","IgG3","IgG4","IgD","IgE")
if (calculate) {
  isotype_usage_timecourse = list()
  for (pat in sclero_pats) {
    file = paste(path,"/",pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab",sep="")
    data0 = readChangeoDb(file)
    data0$SUBISOTYPE = data0$PRIMER
    isotype_usage_timecourse[[pat]] = list()
    for (iso in isotypes) {
      data1 = subset(data0,SUBISOTYPE==iso)
      for (week in sort(unique(data0$WEEK_CORRECTED))) {
        data = subset(data1,WEEK_CORRECTED==week)
        isotype_usage_timecourse[[pat]][[iso]][as.character(week)] = nrow(data)
      }
    }
    saveRDS(isotype_usage_timecourse,paste(path,"subisotype_usage_timecourse.RDS",sep=""))
  }
}

##### Plots
### subisotype usage time courses
isotype_usage_timecourse = readRDS(paste(path,"subisotype_usage_timecourse.RDS",sep=""))
isotypes = c("IgD","IgM","IgA1","IgA2","IgG1","IgG2","IgG3","IgG4","IgE")
iso_colors = brewer.pal(length(isotypes),"Paired")
names(iso_colors) = as.character(isotypes)
pats=armA_pats
for (pat in pats) {
  if (save_eps) {
    pdf(file=paste(fig_path,"/",id_to_dummy[pat],"_subisotype_usage_timecourses.eps",sep=""), 
        width=3.2,height=2,paper="a4",
        colormodel="rgb",pointsize=10)
  }
  par(mar=c(5,6,2,2))
  y=isotype_usage_timecourse[[pat]]
  max_y = max(sapply(y,function(s) if (s[1]==0) {NA} else {max(s)/s[1]}),na.rm=T)
  xlabel="Time (weeks)"
  ylabel="Number of sequences\nrelative to baseline"
  xlimit=c(-1,110)
  ylimit=c(0,8)
  new_plot=T
  for (i in 1:length(names(y))) {
    iso = names(y)[i]
    X = as.numeric(names(y[[iso]]))
    Y = as.vector(y[[iso]])
      Y = Y/Y[1] # normalize to baseline value
      if (new_plot) {
        plot(X[order(X)],Y[order(X)],xlim=xlimit,
             xlab=xlabel,ylab=ylabel,
             ylim=ylimit,log=logaxis,
             axes=F,ann=F,type='o',
             pch=20,col=iso_colors[iso])
        new_plot=F
      } else {
        lines(X[order(X)],Y[order(X)],type='o',
              pch=20,col=iso_colors[iso])
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
  pdf(file=paste(fig_path,"/subisotype_usage_timecourses_legend.eps",sep=""), 
      width=10,height=3,paper="a4",
      colormodel="rgb",pointsize=10)
}
par(mar=c(0,0,0,0))
plot.new()
legend("center",names(iso_colors),pch=rep(15,length(iso_colors)),
       col=iso_colors,bty='n',ncol=5)
if (save_eps) dev.off()