library(data.table)
library(alakazam)
library(RColorBrewer)
extract_Vgene = function(alleles) {
  v=as.vector(alleles)
  w=sapply(v,function(x) {strsplit(x,"\\*")[[1]][1]})
  s=sapply(w,function(x) {paste(strsplit(x,"-")[[1]][1:2],collapse="-")})
  return(as.vector(s))
}

calculate = T
save_eps = T

fig_path = ####
path = ####
sclero_pats = ####
pats = sclero_pats

if (calculate) {
  site_freq_spectrum = list()
  lineage_field = "LINEAGE_dissim0.1"
  lower_size = 5
  upper_size = 30
  fixed_thresh = 0.8
  
  types = c("VERY_SIMILAR","SIMILAR","DISSIMILAR","VERY_DISSIMILAR")
  isotypes = c("IgM","IgD","IgG","IgA")
  mut_field = "NMUT_GERMLINE" 
  
  fmut_fixed_timecourse = percent_radicalmut_timecourse = 
    nmut_IgD_timecourse = frac_seq_IgD_timecourse = mean_n_molec_per_seq_timecourse = numeric()
  fgene_timecourse = list()
  
  for (pat in pats) {  
    site_freq_spectrum[[pat]] = list()
    file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-change-table.tab",sep="")
    data0 = readChangeoDb(file)
    weeks = unique(data0$WEEK_CORRECTED)
    if (pat==###SR3) weeks = setdiff(weeks,weeks[which.min(abs(round(weeks)-24))]) #SR3 had only one sequence at week 24 which would break the script
    for (week in weeks) {
      data = subset(data0,WEEK_CORRECTED==week)
      
      
      ### fixed alleles
      setnames(data,lineage_field,"LINEAGE")
      lineage_sizes = aggregate(COUNT ~ LINEAGE,cbind(data,COUNT=1),FUN=sum)
      lineages = subset(lineage_sizes,COUNT >= lower_size & COUNT < upper_size,LINEAGE,drop=T)
      if (length(lineages)==0) {
        fmut_fixed_timecourse[paste(pat,week,sep="_")] = NA
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
          fmut_fixed_timecourse[paste(pat,week,sep="_")] = 0
        } else {
          fmut_fixed_timecourse[paste(pat,week,sep="_")] = mean(fixed_alleles_by_lineage)
        }
      }
      
      
      ### IGHV2-5 gene
      data$Vgene = extract_Vgene(data$V_CALL_GENOTYPED)
      for (gene in c("IGHV2-5")) {
        ## lineage-weighted:
        linwdata = data[order(data$LINEAGE),]
        linwdata = linwdata[!duplicated(linwdata$LINEAGE),]
        fgene_timecourse[[gene]][paste(pat,week,sep="_")] = sum(linwdata$Vgene==gene)/nrow(linwdata)
      }
        
      
      ### 
      data$ISOTYPE = gsub("1|2|3|4","",data$PRIMER)
      nmut_IgD_timecourse[paste(pat,week,sep="_")] = mean(data[data$ISOTYPE=="IgD",mut_field])
      frac_seq_IgD_timecourse[paste(pat,week,sep="_")] = sum(grepl("IgD",data$PRIMER))/nrow(data)
    }
  }
  
  sclero_sig_timecourse = data.frame(IGHV2.5=fgene_timecourse[["IGHV2-5"]],
                                     fmut_fixed=fmut_fixed_timecourse,
                                     nmut_IgD=nmut_IgD_timecourse,
                                     frac_seq_IgD=frac_seq_IgD_timecourse)
saveRDS(sclero_sig_timecourse,paste(path,"sclero_sig_timecourse.RDS",sep=""))
}

# Make Fig 3 panels
sclero_sig_timecourse = readRDS(paste(path,"sclero_sig_timecourse.RDS",sep=""))
sclero_sig_in_healthy = readRDS(paste(path,"sclero_sig_end-vs-healthy.RDS",sep=""))

colors = brewer.pal(8,"Paired")[c(6,2)]
armA_label = "study arm A"
armB_label = "study arm B"
names(colors) = c(armA_label,armB_label)

armA_pats = ####
features = intersect(colnames(sclero_sig_timecourse),colnames(sclero_sig_in_healthy))
df = rbind(sclero_sig_timecourse[,features],sclero_sig_in_healthy[grepl("healthyBL",rownames(sclero_sig_in_healthy)),features])
df = na.omit(df)
df$IGHV2.5 = 100*df$IGHV2.5
df$fmut_fixed = 100*df$fmut_fixed
df$frac_seq_IgD = 100*df$frac_seq_IgD

for (feature in names(df)) {
  if (save_eps) {
    pdf(file=paste(fig_path,feature,"_timecourse.eps",sep=""), 
        width=3,height=3,paper="a4",
        colormodel="rgb",pointsize = 10)
  }
  par(mar=c(6,6,2,2))
  if (feature=="IGHV2.5") ylabel="% IGHV2-5 usage"
  if (feature=="fmut_fixed") ylabel="% fixed alleles"
  if (feature=="nmut_IgD") ylabel="Mean number of mutations (IgD)"
  if (feature=="frac_seq_IgD") ylabel="% IgD"
yrange = range(na.omit(df[,feature]))
if (feature=="fmut_fixed") yrange=c(0,110)
if (feature=="frac_seq_IgD") yrange[2]=17
for (i in 1:length(sclero_pats)) {
  pat = sclero_pats[i]
  samples = grep(pat,rownames(df),value=T)
  weeks = sapply(samples,function(s) as.numeric(strsplit(s,"_")[[1]][2]))
  ord = order(weeks)
  color = ifelse(pat %in% armA_pats,colors[armA_label],colors[armB_label])
  if (i==1) {
    plot(weeks[ord],df[samples[ord],feature],col=color,type='o',pch=20,
         xlim=c(0,85),ylim=yrange,
         xlab = "Time (weeks)",ylab = ylabel,ann=F,axes=F)
    axis(side = 1, at = axTicks(1), label = axTicks(1))
    axis(side = 2, at = axTicks(2), label = axTicks(2),las=1)
    box(which = "plot", lty = "solid")
    mtext(side=1,line=2,"Time (weeks)")
    mtext(side=2,line=2.25,ylabel)
  } else {
    lines(weeks[ord],df[samples[ord],feature],col=color,type='o',pch=20)
  }
}
min_healthy = min(df[grepl('healthy',rownames(df)),feature])
max_healthy = max(df[grepl('healthy',rownames(df)),feature])
abline(h=min_healthy,lty=2)
abline(h=max_healthy,lty=2)
if (save_eps) dev.off()
}

if (save_eps) {
  pdf(file=paste(fig_path,"sclero-sig-timecourse_legend.eps",sep=""), 
      width=4,height=3,paper="a4",
      colormodel="rgb",pointsize = 10)
}
plot.new()
legend("center",names(colors[c(armA_label,armB_label)]),
       text.col=colors[c(armA_label,armB_label)],bty='n')
if (save_eps) dev.off()