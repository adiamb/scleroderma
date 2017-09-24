library(data.table)
library(alakazam)
library(plyr)
library(tidyr)

calculate = T
save_eps = T

path = ####
fig_path = ####
healthy_pats = ####
sclero_pats = ####

pats = c(healthy_pats,sclero_pats)
lineage_field = "LINEAGE_dissim0.1"
bin_size = 5
fixed_thresh = 0.8

if (calculate) {
  fixed_alleles = list()
  for (i in 1:length(pats)) {  
    pat = pats[i]
    file = paste(path,grep(pat,list.files(path,"attach-change-table.tab"),value=T),sep="")
    data = readChangeoDb(file)
    if ("WEEK" %in% names(data)) data=subset(data,WEEK==0)
    setnames(data,lineage_field,"LINEAGE")
    lineage_sizes = aggregate(COUNT ~ LINEAGE,cbind(data,COUNT=1),FUN=sum)
    bin_breaks = seq(bin_size,max(lineage_sizes$COUNT)+bin_size,bin_size)
    for (i in 1:(length(bin_breaks)-1)) {
      lower_size = bin_breaks[i]
      upper_size = bin_breaks[i+1]
      lineages = subset(lineage_sizes,COUNT >= lower_size & COUNT < upper_size,LINEAGE,drop=T)
      if (length(lineages)==0) {
        fixed_alleles[[pat]][i] = NA
      } else {
        fixed_alleles_by_lineage = sapply(lineages,function(lin) {
          subdata = unique(subset(data,LINEAGE==lin,c("NT_CHANGES","SEQUENCE_IMGT")))
          mutations = unlist(strsplit(subdata$NT_CHANGES,","))
          mutations = mutations[!is.na(mutations)]
          frequencies = table(mutations)/nrow(subdata)
          return(sum(frequencies >= fixed_thresh)/length(frequencies))
        })
        fixed_alleles_by_lineage = fixed_alleles_by_lineage[!is.na(fixed_alleles_by_lineage)]
        if (length(fixed_alleles_by_lineage)==0) {
          fixed_alleles[[pat]][i] = 0
        } else {
          fixed_alleles[[pat]][i] = mean(fixed_alleles_by_lineage)
        }
      }
    }
    saveRDS(fixed_alleles,paste(path,"fixed_alleles_vs_lineage_size.RDS",sep=""))
  }
}


### Fig 2E
if (save_eps) {
  pdf(file=paste(fig_path,"/fixed-alleles-vs-lineage-size_violin-plots.eps",sep=""), 
      width=3.3,height=1.82,paper="a4",
      colormodel="rgb",pointsize=10)
}
set.seed(1)
fixed_alleles = readRDS(paste(path,"fixed_alleles_vs_lineage_size.RDS",sep=""))
group_colors = brewer.pal(6,"Paired")[c(2,6)]
names(group_colors) = c("healthy","scleroderma")
nbins=5
bin_names=paste(bin_size*(1:nbins),"-",bin_size*(2:(nbins+1))-1,sep="")
L=lapply(fixed_alleles,function(s) {result=s[1:nbins]; names(result)=bin_names; return(result)})
df <- data.frame(matrix(unlist(L),nrow=length(L),byrow=T,dimnames=list(names(L),names(L[[1]]))),check.names=F)
df$patient=rownames(df)
df=gather(df,key="LinSizeBin",value="FixedAlleleFrac",1:nbins)
df=df[order(df$patient),]
df$cohort=ifelse(df$patient %in% sclero_pats,"scleroderma","healthy")
df=na.omit(df)
pval.df=data.frame(LinSizeBin=bin_names,
                   cohort="healthy", # nonsensical, only in to keep ggplot happy, doesn't affect anything
                   p=sapply(bin_names,function(bin) {
                     sclero_values = sapply(sclero_pats,function(pat) L[[pat]][bin])
                     healthy_values = sapply(healthy_pats,function(pat) L[[pat]][bin])
                     pval=wilcox.test(sclero_values,healthy_values)$p.value
                     return(pval)}),
                   n=sapply(bin_names,function(bin) {
                     sclero_values = sapply(sclero_pats,function(pat) L[[pat]][bin])
                     healthy_values = sapply(healthy_pats,function(pat) L[[pat]][bin])
                     return(sum(!is.na(sclero_values))+sum(!is.na(healthy_values)))}),
                   pVal=sapply(bin_names,function(bin) {
                     sclero_values = sapply(sclero_pats,function(pat) L[[pat]][bin])
                     healthy_values = sapply(healthy_pats,function(pat) L[[pat]][bin])
                     pval=wilcox.test(sclero_values,healthy_values)$p.value
                     if (pval >= 0.05) symbol=""
                     if (pval < 0.05 & pval >= 0.01) symbol="*  "
                     if (pval < 0.01 & pval >= 0.001) symbol="**"
                     if (pval < 0.001 & pval >= 0.0001) symbol="***"
                     if (pval < 0.0001 & pval >= 0.00001) symbol="****"
                     return(symbol)}))
p <- ggplot(df,aes(factor(LinSizeBin,levels=bin_names),100*FixedAlleleFrac,fill=factor(cohort),colour=factor(cohort))) +
  labs(x="Lineage size",y="% fixed alleles") +
  guides(fill=FALSE,colour=FALSE) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_blank(),
        panel.border=element_rect(colour="black",fill=NA,size=0.7),
        axis.text=element_text(colour='black',family="Helvetica",size=10),
        axis.title=element_text(colour='black',family="Helvetica",size=10)) +
  geom_violin(aes(fill=factor(cohort)),position=position_dodge(0.5)) + 
  scale_colour_manual(values=group_colors) + 
  scale_fill_manual(values=group_colors) +
  geom_point(position=position_jitterdodge(dodge.width=0.5),color='black',pch=21) +
  geom_text(data=pval.df,aes(y=97,label=pval.df$pVal),col='black',size=5)
print(p)
if (save_eps) dev.off()
# Look at identity of outliers
df$DUMMY=id_to_dummy[sapply(strsplit(df$patient,"_"),function(s) s[1])]
df[order(df$LinSizeBin,-df$FixedAlleleFrac),]