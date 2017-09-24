library(RColorBrewer)

save_eps=T

group_colors = brewer.pal(6,"Paired")[c(2,4,6)]
names(group_colors) = c("Study arm B","healthy","Study arm A") 

armA_pats = ###############
armB_pats = #################
healthy_pats = #############
id_to_dummy =  c(paste("SR",1:length(armA_pats),sep=""),
                 paste("SP",1:length(armB_pats),sep=""),
                 paste("H",1:length(healthy_pats),sep=""))
names(id_to_dummy) = c(armA_pats,armB_pats,healthy_pats)


### Data
fig_path = #####
path = #####
D = read.table(paste(path,"RUBL.Nseq-10000.GERMLINE_IMGT.tab",sep=""))

# Fig 1F
if (save_eps) {
  pdf(file=paste(fig_path,"UniFrac_Bl_End.eps",sep=""), 
      width=1.5,height=2.4,paper="a4",
      colormodel="rgb",pointsize = 10)
}
endweek = 36
armA_D = diag(as.matrix(D[paste(armA_pats,"0",sep="_"),paste(armA_pats,endweek,sep="_")]))
names(armA_D) = armA_pats
armB_D = diag(as.matrix(D[paste(armB_pats,"0",sep="_"),paste(armB_pats,endweek,sep="_")]))
names(armB_D) = armB_pats
data_list = list(armA_D,armB_D)
par(mar=c(2,5,8,3))


df = data.frame(UniFrac=c(armA_D,armB_D),arm=c(rep("A",length(armA_D)),rep("B",length(armB_D))),
                row.names=id_to_dummy[names(c(armA_D,armB_D))])
violin_colors = brewer.pal(6,"Paired")[c(2,6)]
names(violin_colors) = c("B","A")
p <- ggplot(df,aes(factor(arm,levels=c("A","B")),UniFrac,fill=factor(arm),colour=factor(arm))) +
  labs(x="Study arm",y="UniFrac(week 0, week 36)") +
  guides(fill=FALSE,colour=FALSE) +
  geom_violin(aes(fill=factor(arm)),position=position_dodge(0.5)) +
  scale_colour_manual(values=violin_colors[df$arm]) + 
  scale_fill_manual(values=violin_colors[df$arm]) +
  geom_point(position=position_jitterdodge(dodge.width=0.5),color='black',pch=21) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background=element_blank(),axis.line=element_blank(),
        panel.border=element_rect(colour="black",fill=NA,size=0.7),
        axis.text=element_text(colour='black',family="Helvetica",size=11),
        axis.title=element_text(colour='black',family="Helvetica",size=11),
        aspect.ratio=2/1)
print(p)
if (save_eps) dev.off()