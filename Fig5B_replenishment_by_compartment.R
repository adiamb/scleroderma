library(alakazam)
library(plyr)
library(data.table)

armA_pats = ###############
armB_pats = #################
healthy_pats = #############
id_to_dummy =  c(paste("SR",1:length(armA_pats),sep=""),
                 paste("SP",1:length(armB_pats),sep=""),
                 paste("H",1:length(healthy_pats),sep=""))
names(id_to_dummy) = c(armA_pats,armB_pats,healthy_pats)

dummy_to_id =  c(armA_pats,armB_pats,healthy_pats)
names(dummy_to_id) = c(paste("SR",1:length(armA_pats),sep=""),
                 paste("SP",1:length(armB_pats),sep=""),
                 paste("H",1:length(healthy_pats),sep=""))

path = ####
fig_path = ####
save_eps = T
calculate = T
replenishment_by_compartment = list()

if (calculate) {
for (pat in armA_pats) {
  file = paste(path,"/",pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab",sep="")
  
  data = readChangeoDb(file)
  data$INCIDENCE = 1
  data$ISOTYPE = gsub("1|2|3|4","",data$PRIMER)
  data$NAIVE = (data$ISOTYPE %in% c("IgD","IgM") & data$NMUT_GERMLINE == 0)
  data$ANTIGEN_EXPERIENCED = !(data$ISOTYPE %in% c("IgD","IgM") & data$NMUT_GERMLINE == 0)
  
  replenishment_by_compartment[[pat]] = aggregate(cbind(INCIDENCE,DUPCOUNT) ~ WEEK_CORRECTED+NAIVE+ANTIGEN_EXPERIENCED, data=data, FUN=sum) 
  
}
saveRDS(replenishment_by_compartment,paste(path,"replenishment_by_compartment.RDS",sep=""))
}

########
replenishment_by_compartment = readRDS(paste(path,"replenishment_by_compartment.RDS",sep=""))
compartment_colors = c(brewer.pal(6,"Paired")[c(2,6)],"grey50")
compartment_lty = c(1,1,1) 
names(compartment_colors) = names(compartment_lty) = c("naive","antigen-experienced","total")

for (pat in armA_pats) {
  data = replenishment_by_compartment[[pat]]
  
  naive = subset(data,NAIVE==T)
  naive = naive[order(naive$WEEK_CORRECTED),]
  antigen_experienced = subset(data,ANTIGEN_EXPERIENCED==T)
  antigen_experienced = antigen_experienced[order(antigen_experienced$WEEK_CORRECTED),]
  total = aggregate(cbind(INCIDENCE,DUPCOUNT) ~ WEEK_CORRECTED,data,FUN=sum)
  total = total[order(total$WEEK_CORRECTED),]
  x_naive = naive$WEEK_CORRECTED
  x_agexp = antigen_experienced$WEEK_CORRECTED
  x_tot = total$WEEK_CORRECTED
  y_naive = naive$DUPCOUNT/total$DUPCOUNT[1]
  y_agexp = antigen_experienced$DUPCOUNT/total$DUPCOUNT[1]
  y_tot = total$DUPCOUNT/total$DUPCOUNT[1]
  if (save_eps) {
    pdf(file=paste(fig_path,id_to_dummy[pat],"_replenishment_by_compartment_v2.eps",sep=""), 
        width=3.2,height=2.3,paper="a4",
        colormodel="rgb",pointsize = 10)
  }
  ylimit=c(0,1.5)
  if (pat %in% ###) ylimit=c(0,4.5)
  par(mar=c(5,6,2,2))
  if (id_to_dummy[pat] %in% c("SR4","SR6")) {xtick_labels=T; xlabel="Time (weeks)"} else {xtick_labels=F; xlabel=""}
  plot(x_tot,y_tot,type='o',
       col=compartment_colors["total"],pch=20,lty=compartment_lty["total"],
       ylim=ylimit,
       xlim=c(-1,110),ann=F,axes=F)
  lines(x_agexp,y_agexp,type='o',pch=20,
       col=compartment_colors["antigen-experienced"],lty=compartment_lty["antigen-experienced"])
  lines(x_naive,y_naive,type='o',pch=20,lty=compartment_lty["naive"],
        col=compartment_colors["naive"])
  if (id_to_dummy[pat] %in% c("SR4","SR6","SR3")) {xtick_labels=T; xlabel="Time (weeks)"} else {xtick_labels=F; xlabel=""}
  mtext(side=1,line=2,xlabel)
  axis(side=1,at=axTicks(1),labels=xtick_labels)
  axis(side=2,at=axTicks(2),las=2,
       label=axTicks(2)
       )
  box()
  abline(1,0,lty=2)
  legend_pos = "topright"
  legend(legend_pos,id_to_dummy[pat],bty='n')
  if (save_eps) dev.off()
}

if (save_eps) {
  pdf(file=paste(fig_path,"replenishment_by_compartment_legend_v2.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend(0,0.5,bty='n',c('naive','antigen-experienced','total'),
       text.col=compartment_colors[c('naive','antigen-experienced','total')])
if (save_eps) dev.off()