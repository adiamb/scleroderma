library(RColorBrewer)

path = ##############
fig_path = ############
save_eps = T
calculate = T


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

  

if (calculate) {
nseq_ratio = fseqswitched_ratio = fsequnswitched_ratio = unmutfseqinunswitched_ratio = numeric()
depl_weeks = c(4,12)
for (pat in sclero_pats) {
  file = paste(path,"/",pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab",sep="")
  data = readChangeoDb(file)
  data$ISOTYPE = gsub("1|2|3|4","",data$PRIMER)
  data$INCIDENCE = 1
  data$SWITCHED_INCIDENCE = as.numeric(data$ISOTYPE %in% c("IgG","IgA"))
  data$UNSWITCHED_INCIDENCE = as.numeric(data$ISOTYPE %in% c("IgM"))
  data$UNMUT_UNSWITCHED_INCIDENCE = as.numeric(data$ISOTYPE %in% c("IgD","IgM") & data$NMUT_GERMLINE==0)
  
  sum_data = aggregate(cbind(INCIDENCE,SWITCHED_INCIDENCE,UNSWITCHED_INCIDENCE,UNMUT_UNSWITCHED_INCIDENCE) ~ WEEK,data=data,FUN=sum)
  sum_data$SWITCHED_FRACTION = sum_data$SWITCHED_INCIDENCE/sum_data$INCIDENCE
  sum_data$UNSWITCHED_FRACTION = sum_data$UNSWITCHED_INCIDENCE/sum_data$INCIDENCE
  sum_data$UNMUT_FRAC_IN_UNSWITCHED = sum_data$UNMUT_UNSWITCHED_INCIDENCE/sum_data$UNSWITCHED_INCIDENCE
  
  nseq_ratio[pat] = mean(subset(sum_data,WEEK %in% depl_weeks,INCIDENCE,drop=T))/subset(sum_data,WEEK==0,INCIDENCE,drop=T)
  fsequnswitched_ratio[pat] = mean(subset(sum_data,WEEK %in% depl_weeks,UNSWITCHED_FRACTION,drop=T))/subset(sum_data,WEEK==0,UNSWITCHED_FRACTION,drop=T)
  fseqswitched_ratio[pat] = mean(subset(sum_data,WEEK %in% depl_weeks,SWITCHED_FRACTION,drop=T))/subset(sum_data,WEEK==0,SWITCHED_FRACTION,drop=T)
  unmutfseqinunswitched_ratio[pat] = mean(subset(sum_data,WEEK %in% depl_weeks,UNMUT_FRAC_IN_UNSWITCHED,drop=T))/subset(sum_data,WEEK==0,UNMUT_FRAC_IN_UNSWITCHED,drop=T)
}
saveRDS(nseq_ratio,paste(path,"nseq_ratio.RDS"))
saveRDS(fseqswitched_ratio,paste(path,"fseqswitched_ratio.RDS"))
saveRDS(fsequnswitched_ratio,paste(path,"fsequnswitched_ratio.RDS"))
saveRDS(unmutfseqinunswitched_ratio,paste(path,"unmutfseqinunswitched_ratio.RDS"))
}

####
nseq_ratio = readRDS(paste(path,"nseq_ratio.RDS"))
fseqswitched_ratio = readRDS(paste(path,"fseqswitched_ratio.RDS"))
fsequnswitched_ratio = readRDS(paste(path,"fsequnswitched_ratio.RDS"))
unmutfseqinunswitched_ratio = readRDS(paste(path,"unmutfseqinunswitched_ratio.RDS"))

### Fig 1E
if (save_eps) {
  pdf(file=paste(fig_path,"depletion_signature.eps",sep=""), 
      width=2.75,height=3,paper="a4",
      colormodel="rgb",pointsize=10)
}
par(mar=c(8,6,2,4))
plot(nseq_ratio[sclero_pats],fseqswitched_ratio[sclero_pats],
     col = "black", bg = sapply(sclero_pats,pat_to_group_color), pch=21,
     log='xy',xlim=c(1e-3,5),ylim=c(0.5,50),
     xlab="",ylab="",xaxt="n",yaxt="n")
at_x <- c(1e-3,1e-2,1e-1,1,10)
at_y <- c(1,2,5,10,20,50)
axis(side=1,at=at_x,las=1,labels=as.expression(sciNotation(at_x)))
axis(side=2,at=at_y,las=2)
mtext(text="Fold-change, sequence diversity",side=1,line=2.1)
mtext(text="Fold-change,\nisotype-switched proportion",side=2,line=2.1)
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"depletion_signature_legend.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend(0.5,0.5,bty='n',c('Study arm A','Study arm B'),
       text.col=group_colors[c('Study arm A','Study arm B')])
if (save_eps) dev.off()