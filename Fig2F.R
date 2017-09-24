tabulate = T
analyze = T
save_eps = T

### Functions
library(alakazam)
library(DESeq2)
library(pheatmap)
extract_Vgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  s=sapply(w,function(x) {paste(strsplit(x,"-")[[1]][1:2],collapse="-")})
  return(as.vector(s))
}
extract_Jgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  w=gsub("Homsap ","",w)
  w[grepl("IGHJ",w)]=NA
  return(as.vector(w))
}
extract_Dgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  w=gsub("Homsap ","",w)
  w[grepl("IGHD",w)]=NA
  return(as.vector(w))
}
correct_Vnames = function(freqs) {
  ### make sure to have consistent gene/allele labeling
  require(ape)
  a = character()
  for (pat in pats) {
    ref1=paste(path,pat,"_personal_hIGHV.fasta",sep="")
    a1=ape::read.FASTA(ref1)
    a1=as.character(a1)
    a1=sapply(a1,function(s) paste(s,collapse=""))
    a = c(a,a1)
  }
  a=gsub("-",".",a)
  a=sapply(1:length(a),function(s) {
    allele = names(a)[s]
    sequence = a[s]
    if (grepl("IGHV1",allele)) {start=277}
    if (grepl("IGHV2",allele)) {start=280}
    if (grepl("IGHV3",allele)) {start=284}
    if (grepl("IGHV4",allele)) {start=281}
    if (grepl("IGHV5",allele)) {start=279}
    if (grepl("IGHV6",allele)) {start=281}
    if (grepl("IGHV7",allele)) {start=288}
    return(substr(sequence,start,nchar(sequence)))  
  }) 
  map = character()
  for (i in 1:length(a)) {
    allele = names(a)[i]
    seq = a[allele]
    map[allele] = sort(names(a)[a==seq])[1]
  }
  names_map = as.character(extract_Vgene(names(map)))
  map = extract_Vgene(map)
  names(map) = names_map
  map2 = character()
  for (name in unique(names(map))) {
    map2[name] = sort(map[names(map)==name])[1]
  }
  for (i in 1:length(map2)) {
    map2[i] = map2[map2[i]]
  }
  colnames(freqs) = map2[colnames(freqs)]
  x = matrix(nrow=nrow(freqs),ncol=length(unique(colnames(freqs))),
             dimnames=list(rownames(freqs),unique(colnames(freqs))))
  for (gene in unique(colnames(freqs))) {
    x[,gene] = rowSums(freqs[,colnames(freqs)==gene,drop=F])
  }
  return(x)
}

### Participants
path = ####
healthy_pats = ####
sclero_pats = ####
pats = c(healthy_pats,sclero_pats)

if (tabulate) {
### Data
all_data = data.frame()
isotypes = c("IgD","IgM","IgA","IgG","all")
for (pat in pats) {
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab",sep="")
  data0 = readChangeoDb(file)
  for (iso in isotypes) {
    if (!(iso=="all")) data = subset(data0,grepl(iso,PRIMER)) else data=data0
    
    data$INCIDENCE = 1
    data$SAMPLE = paste(pat,iso,sep=";")
    data$Vgene = extract_Vgene(data$V_CALL_GENOTYPED)
    
    ### to get lineage-weighted counts:
    data = data[order(data$LINEAGE_dissim0.1),]
    data$COUNT_lindissim0.1 = 1
    data$COUNT_lindissim0.1[duplicated(data$LINEAGE_dissim0.1)] = 0
    ###
    
    all_data = rbind(all_data,data[,c("COUNT_lindissim0.1","SAMPLE","Vgene")])
    
  }
}
Vfreqs_lineages0.1 = xtabs(COUNT_lindissim0.1 ~ SAMPLE + Vgene, all_data)
saveRDS(Vfreqs_lineages0.1,paste(path,"Vfreqs_lineages0.1.RDS",sep=""))
}


if (analyze) {
### DESeq analysis
iso = "all"
freqs = as.data.frame.matrix(readRDS(paste(path,"Vfreqs_lineages0.1.RDS",sep="")))
freqs = freqs[grepl(iso,rownames(freqs)),]
rownames(freqs) = sapply(rownames(freqs), function(s) strsplit(s,";")[[1]][1])
colData = data.frame(SCLERO=as.character(pats %in% sclero_pats),row.names=pats)
dds = DESeqDataSetFromMatrix(countData = t(freqs[pats,]),
                             colData = colData[pats,,drop=F],
                             design = ~ SCLERO)
dds = DESeq(dds)
res = results(dds,contrast=c("SCLERO","TRUE","FALSE"))
res = res[order(res$padj),]
res
diff_genes_deseq = rownames(subset(res,padj<0.15))
}


### Fig 2F
path = ####
healthy_pats = ####
sclero_pats = ####
pats = c(healthy_pats,sclero_pats)
gene = "IGHV2-5"
  if (save_eps) {
    pdf(file=paste(fig_path,gene,"_bySwitchCompartment.eps",sep=""), 
        width=1.75,height=2,paper="a4",
        colormodel="rgb",pointsize=10)
  }
  set.seed(1)
  group_colors = brewer.pal(6,"Paired")[c(2,6)]
  names(group_colors) = c("healthy","scleroderma")
  
  freqs = as.data.frame.matrix(readRDS(paste(path,"Vfreqs_lineages0.1.RDS",sep="")))
  relfreqs = apply(freqs,1,function(s) s/sum(s))
  df = data.frame()
  for (pat in pats) {
      df[pat,"unswitched"] = relfreqs[gene,paste(pat,"IgD",sep=";")]+relfreqs[gene,paste(pat,"IgM",sep=";")]
      df[pat,"switched"] = relfreqs[gene,paste(pat,"IgA",sep=";")]+relfreqs[gene,paste(pat,"IgG",sep=";")]
  }
  df=gather(df,key="compartment",value="geneUsage")
  df$patient=pats
  df$cohort=ifelse(df$patient %in% sclero_pats,"scleroderma","healthy")
  compartments = c("unswitched","switched")
  pval.df=data.frame(compartment=compartments,
                     cohort="healthy", # nonsensical, only in to keep ggplot happy, doesn't affect anything
                     p=sapply(compartments,function(comp) {
                       sclero_values = subset(df,compartment==comp & patient %in% sclero_pats)$geneUsage
                       healthy_values = subset(df,compartment==comp & patient %in% healthy_pats)$geneUsage
                       return(wilcox.test(sclero_values,healthy_values)$p.value)}),
                     n=sapply(compartments,function(comp) {
                       sclero_values = subset(df,compartment==comp & patient %in% sclero_pats)$geneUsage
                       healthy_values = subset(df,compartment==comp & patient %in% healthy_pats)$geneUsage
                       return(sum(!is.na(sclero_values))+sum(!is.na(healthy_values)))}),
                     pVal=sapply(compartments,function(comp) {
                       sclero_values = subset(df,compartment==comp & patient %in% sclero_pats)$geneUsage
                       healthy_values = subset(df,compartment==comp & patient %in% healthy_pats)$geneUsage
                       pval=wilcox.test(sclero_values,healthy_values)$p.value
                       if (pval >= 0.05) symbol=""
                       if (pval < 0.05 & pval >= 0.01) symbol="*"
                       if (pval < 0.01 & pval >= 0.001) symbol="   **"
                       if (pval < 0.001 & pval >= 0.0001) symbol="***"
                       if (pval < 0.0001 & pval >= 0.00001) symbol="****"
                       return(symbol)}))
  p <- ggplot(df,aes(factor(compartment,levels=compartments),100*geneUsage,fill=factor(cohort),colour=factor(cohort))) +
    labs(x="",y=paste("% ",gene," usage",sep="")) +
    guides(fill=FALSE,colour=FALSE) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.background=element_blank(),axis.line=element_blank(),#axis.line=element_line(colour="black"),
          panel.border=element_rect(colour="black",fill=NA,size=0.7),
          axis.text=element_text(colour='black',family="Helvetica",size=10),
          axis.title=element_text(colour='black',family="Helvetica",size=10)) +
    geom_violin(aes(fill=factor(cohort)),position=position_dodge(0.75),width=0.75) + scale_y_log10() +
    scale_colour_manual(values=group_colors) + 
    scale_fill_manual(values=group_colors) +
    scale_x_discrete(labels=c("unswitched" = "IgD/\nIgM", "switched" = "IgG/\nIgA")) +
    geom_point(position=position_jitterdodge(dodge.width=0.75),color='black',pch=21) +
    geom_text(data=pval.df,aes(y=max((90)*df$geneUsage),label=pval.df$pVal),col='black',size=5)
  print(p)
  if (save_eps) dev.off()