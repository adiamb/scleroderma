library(RColorBrewer)
library(pheatmap)
library(alakazam)
library(data.table)

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

log_delta = 10
lineage_field = "LINEAGE_dissim0.1"

if (calculate) {
  size_usage_timecourse = list()
  for (pat in sclero_pats) {
    file = paste(path,"/",pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab",sep="")
    data0 = readChangeoDb(file)
    setnames(data0,lineage_field,"LINEAGE")
    data0 = cbind(unique(data0[,c("LINEAGE","SEQUENCE_IMGT","WEEK_CORRECTED")]),COUNT=1)
    data0 = aggregate(COUNT ~ LINEAGE+WEEK_CORRECTED,data0,sum)
    data0$SIZE_BINNED = log_delta^(ceiling(log(data0$COUNT)/log(log_delta)))
    size_usage_timecourse[[pat]] = list()
    for (size in sort(unique(data0$SIZE_BINNED))) {
      data1 = subset(data0,SIZE_BINNED==size)
      for (week in sort(unique(data0$WEEK_CORRECTED))) {
        data = subset(data1,WEEK_CORRECTED==week)
        size_usage_timecourse[[pat]][[as.character(size)]][as.character(week)] = nrow(data)
      }
    }
    saveRDS(size_usage_timecourse,paste(path,"lineage-size_usage_timecourse.RDS",sep=""))
    }
  }

##### Plots
size_colors = brewer.pal(8,"Paired")[c(2,4,6)]
names(size_colors) = c("1","2-10","> 10")

size_usage_timecourse = readRDS(paste(path,"lineage-size_usage_timecourse.RDS",sep=""))
pats = names(size_usage_timecourse)
for (pat in pats) {
  width=3.2
  height=2
  if (save_eps) {
    pdf(file=paste(fig_path,"/",id_to_dummy[pat],"_lineage-size_usage_threshold_timecourses.eps",sep=""), 
        width=width,height=height,paper="a4",
        colormodel="rgb",pointsize=10)
  }
  par(mar=c(5,6,2,2))
  y=size_usage_timecourse[[pat]]
  if (length(y) < 4) y[[4]]=0
  y2 = list(y[[1]],y[[2]],y[[3]]+y[[4]]); names(y2) =  c("1","2 to 10","> 10"); y = y2
  max_y = max(sapply(y,function(s) if (s[1]==0) {NA} else {max(s)/s[1]}),na.rm=T)
  xlabel="Time (weeks)"
  ylabel="Number of lineages\nrelative to baseline"
  xlimit=c(-1,110)
  if (pat==###SR3) ylimit=c(0,9) else ylimit=c(0,3)
  logaxis=''
  new_plot=T
  for (i in 1:length(names(y))) {
    size = names(y)[i]
    X = as.numeric(names(y[[size]]))
    Y = as.vector(y[[size]])
    Y=Y/Y[1]
    if (new_plot) {
      plot(X[order(X)],Y[order(X)],xlim=xlimit,
           xlab=xlabel,ylab=ylabel,
           ylim=ylimit,log=logaxis,
           axes=F,ann=F,type='o',
           pch=20,cex=1,col=size_colors[size])
      new_plot=F
    } else {
      lines(X[order(X)],Y[order(X)],type='o',pch=20,cex=1,col=size_colors[size])
    }
  }
  if (id_to_dummy[pat] %in% c("SR6")) {xtick_labels=T; xlabel="Time (weeks)"} else {xtick_labels=F; xlabel=""}
  axis(side = 1, at = axTicks(1), labels = xtick_labels)
  if (id_to_dummy[pat]=="SR3") {
    axis(side=2,at=axTicks(2),label=axTicks(2),las=1)
  } else {
    axis(side=2,at=0:3,label=0:3,las=1)
  }
  box(which = "plot", lty = "solid")
  mtext(side = 1, text = xlabel, line = 2)
  abline(1,0,lty=2)
  legend_pos="topright"
  if (id_to_dummy[pat]=="SR6") {legend_pos="topleft"}
  legend(legend_pos,id_to_dummy[pat],bty='n')
  if (save_eps) dev.off()
}

if (save_eps) {
  pdf(file=paste(fig_path,"/lineage-size_usage_threshold_timecourses_legend.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize=10)
}
par(mar=c(0,0,0,0))
plot.new()
legend("center",names(size_colors),col=size_colors,pch=rep(15,length(size_colors)),bty='n',ncol=3)
if (save_eps) dev.off()
