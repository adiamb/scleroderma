library(RColorBrewer)

sciNotation <- function(x, digits = 1) {
  if (length(x) > 1) {
    return(append(sciNotation(x[1]), sciNotation(x[-1])))
  }
  if (!x) return(0)
  exponent <- floor(log10(x))
  base <- round(x / 10^exponent, digits)
  if (base==1) {substitute(10^exponent,
                           list(base = base, exponent = exponent))}
  else {substitute(base %*% 10^exponent,
                   list(base = base, exponent = exponent))}
} 


armA_pats = ###############
armB_pats = #################
healthy_pats = #############
id_to_dummy =  c(paste("SR",1:length(armA_pats),sep=""),
                 paste("SP",1:length(armB_pats),sep=""),
                 paste("H",1:length(healthy_pats),sep=""))
names(id_to_dummy) = c(armA_pats,armB_pats,healthy_pats)


group_colors = brewer.pal(6,"Paired")[c(2,4,6)]
names(group_colors) = c("Study arm B","healthy","Study arm A")
pat_to_group_color = function(s) {
  if (s %in% armA_pats) return(group_colors["Study arm A"])
  if (s %in% healthy_pats) return(group_colors["healthy"])
  if (s %in% armB_pats) return(group_colors["Study arm B"])
}

path = ####
fig_path = ####
save_eps = T
calculate = T

sclero_pats = ####

if (calculate) {
n_seq = n_molec = 
  frac_seq_switched = unmut_frac_seq_in_unswitched = list()
for (pat in sclero_pats) {
  file = paste(path,"/",pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab",sep="")
  data0 = readChangeoDb(file)
  data0$ISOTYPE = gsub("1|2|3|4","",data0$PRIMER)
  for (week in sort(unique(data0$WEEK_CORRECTED))) {
    data = subset(data0,WEEK_CORRECTED==week)
    n_seq[[pat]][as.character(week)] = nrow(data)
    n_molec[[pat]][as.character(week)] = sum(data$DUPCOUNT)
    f=data$DUPCOUNT/sum(data$DUPCOUNT)
    frac_seq_switched[[pat]][as.character(week)] = sum(data$ISOTYPE %in% c("IgG","IgA"))/nrow(data)
    unmut_frac_seq_in_unswitched[[pat]][as.character(week)] = sum(data$NMUT_GERMLINE==0 & data$ISOTYPE %in% c("IgD","IgM"))/sum(data$ISOTYPE %in% c("IgD","IgM"))
  }
}

saveRDS(n_seq,paste(path,"n_seq","_timecourse.RDS"))
saveRDS(n_molec,paste(path,"n_molec","_timecourse.RDS"))
saveRDS(frac_seq_switched,paste(path,"frac_seq_switched","_timecourse.RDS"))
saveRDS(unmut_frac_seq_in_unswitched,paste(path,"unmut_frac_seq_in_unswitched","_timecourse.RDS"))
}

#####
n_seq = readRDS(paste(path,"n_seq","_timecourse.RDS"))
n_molec = readRDS(paste(path,"n_molec","_timecourse.RDS"))
frac_seq_switched = readRDS(paste(path,"frac_seq_switched","_timecourse.RDS"))
unmut_frac_seq_in_unswitched = readRDS(paste(path,"unmut_frac_seq_in_unswitched","_timecourse.RDS"))


# Fig 1B
if (save_eps) {
  pdf(file=paste(fig_path,"nseq_timecourse.eps",sep=""), 
      width=3.2,height=3,paper="a4",
      colormodel="rgb",pointsize=10)
}
par(mar=c(8,6,2,2))
y=n_seq
xlabel="Time (weeks)"
ylabel="Observed diversity"
xlimit=c(-1,110)
ylimit=c(1,max(unlist(y)))
logaxis='y'
cex=1
for (i in 1:length(names(y))) {
  pat = names(y)[i]
  X = as.numeric(names(y[[pat]]))
  Y = as.vector(y[[pat]])
  if (i==1) {
    plot(X[order(X)],Y[order(X)],xlim=xlimit,
         xlab=xlabel,ylab=ylabel,
         ylim=ylimit,log=logaxis,
         axes=F,ann=F,type='o',
         pch=20,cex=cex,col=pat_to_group_color(pat))
  } else {
    lines(X[order(X)],Y[order(X)],type='o',
          pch=20,cex=cex,col=pat_to_group_color(pat))
  } 
}
axis(side = 1, at = axTicks(1), label = axTicks(1))
axis(side = 2, at = axTicks(2), label = as.expression(sciNotation(axTicks(2),1)),las=1)
box(which = "plot", lty = "solid")
mtext(side = 1, text = xlabel, line = 2.1)
mtext(side = 2, text = ylabel, line = 2.5)
if (save_eps) dev.off()


# Fig 1C
if (save_eps) {
  pdf(file=paste(fig_path,"frac_seq_switched_timecourse.eps",sep=""), 
      width=3.2,height=3,paper="a4",
      colormodel="rgb",pointsize=10)
}
par(mar=c(8,6,2,2))
y=frac_seq_switched
xlabel="Time (weeks)"
ylabel="Isotype-switched proportion"
xlimit=c(-1,110)
ylimit=c(0,max(unlist(y)))
logaxis=''
for (i in 1:length(names(y))) {
  pat = names(y)[i]
  X = as.numeric(names(y[[pat]]))
  Y = as.vector(y[[pat]])
  if (i==1) {
    plot(X[order(X)],Y[order(X)],xlim=xlimit,
         xlab=xlabel,ylab=ylabel,
         ylim=ylimit,log=logaxis,
         axes=F,ann=F,type='o',
         pch=20,cex=cex,col=pat_to_group_color(pat))
  } else {
    lines(X[order(X)],Y[order(X)],type='o',
          pch=20,cex=cex,col=pat_to_group_color(pat))
  }
}
axis(side = 1, at = axTicks(1), label = axTicks(1))
axis(side = 2, at = axTicks(2), label = axTicks(2),las=1)
box(which = "plot", lty = "solid")
mtext(side = 1, text = xlabel, line = 2.1)
mtext(side = 2, text = ylabel, line = 2.5)
if (save_eps) dev.off()


# Fig 1D
if (save_eps) {
  pdf(file=paste(fig_path,"unmut_frac_seq_in_unswitched_timecourse.eps",sep=""), 
      width=3.2,height=3,paper="a4",
      colormodel="rgb",pointsize=10)
}
par(mar=c(8,6,2,2))
y=unmut_frac_seq_in_unswitched
xlabel="Time (weeks)"
ylabel="Unmutated proportion of\nunswitched compartment"
xlimit=c(-1,110)
ylimit=c(0,max(unlist(y),na.rm=T))
logaxis=''
for (i in 1:length(names(y))) {
  pat = names(y)[i]
  X = as.numeric(names(y[[pat]]))
  Y = as.vector(y[[pat]])
  cond = (!is.na(Y))
  X = X[cond]
  Y = Y[cond]
  if (i==1) {
    plot(X[order(X)],Y[order(X)],xlim=xlimit,
         xlab=xlabel,ylab=ylabel,
         ylim=ylimit,log=logaxis,
         axes=F,ann=F,type='o',
         pch=20,cex=cex,col=pat_to_group_color(pat))
  } else {
    lines(X[order(X)],Y[order(X)],type='o',
          pch=20,cex=cex,col=pat_to_group_color(pat))
  }
}
axis(side = 1, at = axTicks(1), label = axTicks(1))
axis(side = 2, at = axTicks(2), label = axTicks(2),las=1)
box(which = "plot", lty = "solid")
mtext(side = 1, text = xlabel, line = 2.1)
mtext(side = 2, text = ylabel, line = 2.5)
if (save_eps) dev.off()