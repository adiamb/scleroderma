calculate = T
save_eps = T

fig_path = ####
path = #####
healthy_pats = ##########
sclero_pats = #############
pats = c(healthy_pats,sclero_pats)

if (calculate) {
  mean_n_molec_per_seq_bl = numeric()

for (pat in pats) {  
  file = paste(path,pat,
"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds
-pass_attach-AA-sequences-pass-corrected.tab",sep="")
  
  data = readChangeoDb(file)

  mean_n_molec_per_seq_bl[pat] = sum(data$DUPCOUNT)/nrow(data)
  
  }
saveRDS(mean_n_molec_per_seq_bl,paste(path,"mean_n_molec_per_seq_bl.RDS",sep=""))
}


mean_n_molec_per_seq_bl = readRDS(paste(path,"mean_n_molec_per_seq_bl.RDS",sep=""))


### Fig 2B
mean_n_molec_per_seq_bl = readRDS(paste(path,"mean_n_molec_per_seq_bl.RDS",sep=""))

if (save_eps) {
  pdf(file=paste(fig_path,"mean-transcripts-per-sequence_violin-pot.eps",sep=""), 
      width=1.5,height=2,paper="a4",
      colormodel="rgb",pointsize=10)
}
set.seed(1)
group_colors = brewer.pal(6,"Paired")[c(2,6)]
names(group_colors) = c("healthy","scleroderma")
df = data.frame(mean_n_molec_per_seq_bl=mean_n_molec_per_seq_bl,row.names=names(mean_n_molec_per_seq_bl))
df$patient=rownames(df)
df$cohort=ifelse(df$patient %in% sclero_pats,"scleroderma","healthy")
df$isotype="all isotypes"
pval.df=data.frame(isotype="all isotypes",cohort="healthy", # nonsensical, only in to keep ggplot happy, doesn't affect anything
                   pVal=wilcox.test(df$mean_n_molec_per_seq_bl[df$cohort=="scleroderma"],
                                    df$mean_n_molec_per_seq_bl[df$cohort=="healthy"])$p.value,
                   n=sum(!is.na(df$mean_n_molec_per_seq_bl[df$cohort=="scleroderma"]))+
                     sum(!is.na(df$mean_n_molec_per_seq_bl[df$cohort=="healthy"])))
pval.df$pSymbol = sapply(pval.df$pVal, function(s) {
  if (s < 0.05 & s >= 0.01) symbol="*"
  if (s < 0.01 & s >= 0.001) symbol="**"
  if (s < 0.001 & s >= 0.0001) symbol="***"
  if (s < 0.0001 & s >= 0.00001) symbol="****"
  return(symbol)
})
p <- ggplot(df,aes(factor(isotype),mean_n_molec_per_seq_bl,fill=factor(cohort),colour=factor(cohort))) +
  labs(x="",y="Mean transcripts\nper IGH sequence",size=6) +
  ylim(0,max(1.05*df$mean_n_molec_per_seq_bl)) +
  guides(fill=FALSE,colour=FALSE) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_blank(),
        panel.border=element_rect(colour="black",fill=NA,size=0.7),
        axis.text=element_text(colour='black',family="Helvetica",size=10),
        axis.title=element_text(colour='black',family="Helvetica",size=10)) +
  geom_violin(aes(fill=factor(cohort)),position=position_dodge(1),width=0.5) + 
  scale_colour_manual(values=group_colors) + 
  scale_fill_manual(values=group_colors) +
  geom_point(position=position_jitterdodge(dodge.width=1),color='black',pch=21) +
  geom_text(data=pval.df,aes(y=max(1*df$mean_n_molec_per_seq_bl),label=pval.df$pSymbol),col='black',size=5)
print(p)
if (save_eps) dev.off()
View(pval.df)