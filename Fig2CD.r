library(RColorBrewer)
library(alakazam)
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

calculate = T
save_eps = T

fig_path = #########
path = ############
healthy_pats = ##############
sclero_pats = ###############
pats = c(healthy_pats,sclero_pats)

if (calculate) {
abund_hist_compartment = list()
fseq_highabund_switched_new = fseq_highabund_unswitched_new = numeric()
bin_size = 20
frac_bin_size = 1e-4
abund_thresh = 5e-4

for (i in 1:length(pats)) {  
  pat = pats[i]
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab",sep="")
  data = readChangeoDb(file)
  data = aggregate(DUPCOUNT ~ SEQUENCE_VDJ+PRIMER,data=data,FUN=sum)
  data$ISOTYPE = gsub("1|2|3|4","",data$PRIMER)
  data$FRAC_ABUNDANCE = data$DUPCOUNT/sum(data$DUPCOUNT)
  abund_hist_compartment[[pat]][["unswitched"]] = hist(subset(data,ISOTYPE %in% c("IgD","IgM"),DUPCOUNT,drop=T),
                                                            breaks = seq(-bin_size/2,3*bin_size/2+max(data$DUPCOUNT),bin_size),plot=F)
  abund_hist_compartment[[pat]][["switched"]] = hist(subset(data,ISOTYPE %in% c("IgA","IgG"),DUPCOUNT,drop=T),
                                                          breaks = seq(-bin_size/2,3*bin_size/2+max(data$DUPCOUNT),bin_size),plot=F)
  
  switched_fracdupcounts = subset(data,ISOTYPE %in% c("IgA","IgG"),FRAC_ABUNDANCE,drop=T)
  unswitched_fracdupcounts = subset(data,ISOTYPE %in% c("IgD","IgM"),FRAC_ABUNDANCE,drop=T)
  fseq_highabund_switched_new[pat] = sum(switched_fracdupcounts >= abund_thresh)/length(switched_fracdupcounts)
  fseq_highabund_unswitched_new[pat] = sum(unswitched_fracdupcounts >= abund_thresh)/length(unswitched_fracdupcounts)
}
saveRDS(abund_hist_compartment,paste(path,"abund_hist_compartment.RDS",sep=""))
saveRDS(fseq_highabund_switched_new,paste(path,"fseq_highabund_switched_new.RDS",sep=""))
saveRDS(fseq_highabund_unswitched_new,paste(path,"fseq_highabund_unswitched_new.RDS",sep=""))
}

### Plot
abund_hist_compartment = readRDS(paste(path,"abund_hist_compartment.RDS",sep=""))
fseq_highabund_switched_new = readRDS(paste(path,"fseq_highabund_switched_new.RDS",sep=""))
fseq_highabund_unswitched_new = readRDS(paste(path,"fseq_highabund_unswitched_new.RDS",sep=""))


# Fig 2C
if (save_eps) {
  pdf(file=paste(fig_path,"abundance-distributions.eps",sep=""), 
      width=3,height=2.5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
par(mar=c(4,4,2,2))
data_list = abund_hist_compartment
colors = brewer.pal(8,"Paired")[c(2,4,6,8)]
healthy_switched_label = "healthy, IgG/IgA"
healthy_unswitched_label = "healthy, IgD/IgM"
scleroderma_switched_label = "SSc-PAH, IgG/IgA"
scleroderma_unswitched_label = "SSc-PAH, IgD/IgM"
names(colors) = c(healthy_switched_label,healthy_unswitched_label,
                  scleroderma_switched_label,scleroderma_unswitched_label)
for (i in 1:length(pats)) {  
  pat = pats[i]
  h = data_list[[pat]][["switched"]]
  x = h$mids
  y = h$counts/sum(h$counts)
  color = ifelse(pat %in% sclero_pats,colors[scleroderma_switched_label],colors[healthy_switched_label])
  if (i==1) {
    plot(x,y,col=color,type='o',log='xy',ann=F,axes=F,ylim=c(1e-6,1),pch=20,cex=0.7)
    mtext(side=1,line=2,"Abundance")
    mtext(side=2,line=3,"Frequency")
    axis(side=2,at=c(1e-6,1e-4,1e-2,1),label=c(as.expression(sciNotation(c(1e-6,1e-4,1e-2))),1),las=2)
    axis(side=1,at=axTicks(1))
    box()
  } else {
    lines(x,y,col=color,type='o',pch=20,cex=0.7)
  }
  h = data_list[[pat]][["unswitched"]]
  x = h$mids
  y = h$counts/sum(h$counts)
  color = ifelse(pat %in% sclero_pats,colors[scleroderma_unswitched_label],colors[healthy_unswitched_label])
  lines(x,y,col=color,type='o',pch=20,cex=0.7)
}
if (save_eps) dev.off()

if (save_eps) {
  pdf(file=paste(fig_path,"abundance-distributions-legend.eps",sep=""), 
      width=10,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend("center",c("IgG/IgA:","healthy","SSc-PAH","IgD/IgM:","healthy","SSc-PAH"),
       text.col=c("black",colors[c(healthy_switched_label,scleroderma_switched_label)],
                  "black",colors[c(healthy_unswitched_label,scleroderma_unswitched_label)]),
       bty='n',ncol=2)
if (save_eps) dev.off()


# Fig 2D
library(ggplot2); library(tidyr)
fseq_highabund_unswitched = fseq_highabund_unswitched_new
fseq_highabund_switched = fseq_highabund_switched_new
if (save_eps) {
  pdf(file=paste(fig_path,"/fseq-highabund_violin-plots.eps",sep=""), 
      width=2,height=2,paper="a4",
      colormodel="rgb",pointsize=10)
}
set.seed(1)
group_colors = brewer.pal(6,"Paired")[c(2,6)]
names(group_colors) = c("healthy","scleroderma")
df = data.frame(unswitched=fseq_highabund_unswitched,
                switched=fseq_highabund_switched,
                patient=pats,cohort=ifelse(pats %in% sclero_pats,"scleroderma","healthy"),
                row.names=pats)
df=gather(df,key="compartment",value="fseq_highabund",c(unswitched,switched))
compartments = unique(df$compartment)
pval.df=data.frame(compartment=compartments,
                   cohort="healthy", # nonsensical, only in to keep ggplot happy, doesn't affect anything
                   p=sapply(compartments,function(comp) {
                     sclero_values = subset(df,patient %in% sclero_pats & compartment==comp)$fseq_highabund
                     healthy_values = subset(df,patient %in% healthy_pats & compartment==comp)$fseq_highabund
                     pval=wilcox.test(sclero_values,healthy_values)$p.value
                     return(pval)}),
                   n=sapply(compartments,function(comp) {
                     sclero_values = subset(df,patient %in% sclero_pats & compartment==comp)$fseq_highabund
                     healthy_values = subset(df,patient %in% healthy_pats & compartment==comp)$fseq_highabund
                     return(sum(!is.na(sclero_values))+sum(!is.na(healthy_values)))}),
                   pVal=sapply(compartments,function(comp) {
                     sclero_values = subset(df,patient %in% sclero_pats & compartment==comp)$fseq_highabund
                     healthy_values = subset(df,patient %in% healthy_pats & compartment==comp)$fseq_highabund
                     pval=wilcox.test(sclero_values,healthy_values)$p.value
                     if (pval < 0.05 & pval >= 0.01) symbol="*"
                     if (pval < 0.01 & pval >= 0.001) symbol="**"
                     if (pval < 0.001 & pval >= 0.0001) symbol="***   "
                     if (pval < 0.0001 & pval >= 0.00001) symbol="****"
                     return(symbol)}))
options(scipen=10000)
p <- ggplot(df,aes(factor(compartment,levels=compartments),100*fseq_highabund,fill=factor(cohort),colour=factor(cohort))) +
  labs(x="",y="% highly represented\nsequences") +
  guides(fill=FALSE,colour=FALSE) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_blank(),
        panel.border=element_rect(colour="black",fill=NA,size=0.7),
        axis.text=element_text(colour='black',family="Helvetica",size=10),
        axis.title=element_text(colour='black',family="Helvetica",size=10)) +
  geom_violin(aes(fill=factor(cohort)),position=position_dodge(0.75),width=0.5) + scale_y_log10() +
  scale_x_discrete(labels=c("unswitched" = "IgD/\nIgM", "switched" = "IgG/\nIgA")) +
  scale_colour_manual(values=group_colors) + 
  scale_fill_manual(values=group_colors) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),color='black',pch=21) +
  geom_text(data=pval.df,aes(y=0.8*max(100*df$fseq_highabund),label=pval.df$pVal),col='black',size=5)
print(p)
if (save_eps) dev.off()