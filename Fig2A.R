calculate = T
save_eps = T
path = #######
healthy_pats = ######
sclero_pats = ########
pats = c(healthy_pats,sclero_pats)

if (calculate) {
nmut_IgD_bl = numeric()
pats = c(healthy_pats,sclero_pats)
mut_field = "NMUT_GERMLINE"
for (pat in pats) {  
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab",sep="") # data file for baseline
  data0 = readChangeoDb(file)
  data0$ISOTYPE = gsub("1|2|3|4","",data0$PRIMER)
  nmut_IgD_bl[pat] = mean(data0[data0$ISOTYPE=="IgD",mut_field])
  frac_seq_IgD_bl[pat] = sum(data0$ISOTYPE=="IgD")/nrow(data0)
}
saveRDS(nmut_IgD_bl,paste(path,"nmut_IgD_bl.RDS",sep=""))
saveRDS(frac_seq_IgD_bl,paste(path,"frac_seq_IgD_bl.RDS",sep=""))
}


# Fig 2A
nmut_IgD_bl = readRDS(paste(path,"nmut_IgD_bl.RDS",sep=""))
frac_seq_IgD_bl = readRDS(paste(path,"frac_seq_IgD_bl.RDS",sep=""))
group_colors = brewer.pal(6,"Paired")[c(2,6)]
healthy_label = "healthy"
sclero_label = "SSc-PAH"
names(group_colors) = c(healthy_label,sclero_label) 
if (save_eps) {
  pdf(file=paste(fig_path,"IgD-fracseq-vs-nmut.eps",sep=""), 
      width=2.5,height=2.5,paper="a4",
      colormodel="rgb",pointsize=10)
}
par(mar=c(4,4,2,2))
colors = ifelse(pats %in% sclero_pats,group_colors[sclero_label],group_colors[healthy_label])
plot(100*frac_seq_IgD_bl[pats],nmut_IgD_bl[pats],
     col = "black", bg = colors,
     pch=21,log='xy',
     xlab="",ylab="")
mtext(text="% IgD",side=1,line=2)
mtext(text="Mean number of mutations (IgD)",side=2,line=2)
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"IgD-fracseq-vs-nmut_legend.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend(0.5,0.5,bty='n',c(healthy_label,sclero_label),
       text.col=group_colors[c(healthy_label,sclero_label)])
if (save_eps) dev.off()