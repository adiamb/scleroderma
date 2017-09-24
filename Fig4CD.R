save_eps = T
fig_path = ####

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

replenishment_by_compartment = readRDS(paste(path,"replenishment_by_compartment.RDS",sep="")) # made using script Fig5B_replenishment_by_compartment.R
t1 = t2 = recovery_time = initial_diversity = final_diversity = 
  initial_diversity_all = final_diversity_all = numeric()
for (pat in names(replenishment_by_compartment)) {
  data = replenishment_by_compartment[[pat]]

  alldata = aggregate(cbind(INCIDENCE,DUPCOUNT) ~ WEEK_CORRECTED,data,FUN=sum)
  alldata = alldata[order(alldata$WEEK_CORRECTED),]
  initial_diversity_all[pat] = alldata$INCIDENCE[1]
  final_diversity_all[pat] = subset(alldata,WEEK_CORRECTED > 30 & WEEK_CORRECTED < 40,INCIDENCE,drop=T)

  data = subset(data,NAIVE==T)
  data = data[order(data$WEEK_CORRECTED),]
  initial_diversity[pat] = data$INCIDENCE[1]
  final_diversity[pat] = subset(data,WEEK_CORRECTED > 30 & WEEK_CORRECTED < 40,INCIDENCE,drop=T)
  
  t1[pat] = data$WEEK_CORRECTED[2]
  end_count = data$INCIDENCE[length(data$INCIDENCE)]
  replenished = subset(data[-c(1,2),],INCIDENCE > 5e3,WEEK_CORRECTED,drop=T)
  t2[pat] = replenished[1]
  recovery_time[pat] = t2[pat]-t1[pat]
}


# Fig 4D
if (save_eps) {
  pdf(file=paste(fig_path,"initial-diversity-vs-recovery-time.eps",sep=""), 
      width=3,height=2.5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
x = initial_diversity
y = recovery_time
xlabel = "Naive diversity, week 0"
ylabel = "Time to onset\nof replenishment (weeks)"
par(mar=c(6,6,2,2))
plot(x,y,ann=F,axes=F,
     xlab=xlabel,ylab=ylabel,pch=20)
corrtest = cor.test(x,y)
mtext(side=1,line=2,xlabel)
mtext(side=2,line=2,ylabel)
xticks = seq(0,1.5e5,5e4)
axis(side=1,at=xticks,label=as.expression(sciNotation(xticks, 1)))
axis(side=2,at=axTicks(2),las=2)
box()
legend("topright",paste("correlation: ",round(corrtest$estimate,2),
                       "\np-value: ",round(corrtest$p.value,4),
                       "\n(PPMC, n=",length(x),")",sep=""),
       cex=0.8,bty='n')
if (save_eps) dev.off()


# initial diversity correlated with final diversity
if (save_eps) {
  pdf(file=paste(fig_path,"initial-vs-final-diversity.eps",sep=""), 
      width=3,height=2.5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
xlabel = "Naive diversity, week 0"
ylabel = "Naive diversity, week 36"
par(mar=c(6,6,2,2))
plot(initial_diversity,final_diversity,ann=F,axes=F,
     xlab=xlabel,ylab=ylabel,pch=20)
corrtest = cor.test(initial_diversity,final_diversity)
mtext(side=1,line=2,xlabel)
mtext(side=2,line=4.5,ylabel)
xticks = seq(0,1.5e5,5e4)
axis(side=1,at=xticks,label=as.expression(sciNotation(xticks, 1)))
axis(side=2,at=axTicks(2),las=2,label=as.expression(sciNotation(axTicks(2), 1)))
box()
legend("topleft",paste("correlation = ",round(corrtest$estimate,2),
                       "\np = ",round(corrtest$p.value,4),
                       "\n(PPMC, n=",length(initial_diversity),")",sep=""),
       cex=0.8,bty='n')
if (save_eps) dev.off()


# Fig 4C
fit_a_chao = readRDS(paste(path,"a_constrainedfit","_chao.RDS",sep=""))
fit_b_chao = readRDS(paste(path,"b_constrainedfit","_chao.RDS",sep=""))
fit_c1_chao = readRDS(paste(path,"c1_constrainedfit","_chao.RDS",sep=""))
fit_c2_chao = readRDS(paste(path,"c2_constrainedfit","_chao.RDS",sep=""))

if (save_eps) {
  pdf(file=paste(fig_path,"initial-diversity-vs-constrainedfit-a.eps",sep=""), 
      width=3,height=2.5,paper="a4",
      colormodel="rgb",pointsize = 10)
}
fit_pats = names(fit_a_chao)
x = initial_diversity[fit_pats]
y = fit_a_chao[fit_pats]
xlabel = "Naive diversity, week 0"
ylabel = "a (1/week)"
par(mar=c(6,6,2,2))
plot(x,y,ann=F,axes=F,
     xlab=xlabel,ylab=ylabel,pch=20)
corrtest = cor.test(x,y)
mtext(side=1,line=2,xlabel)
mtext(side=2,line=4.5,ylabel)
xticks = seq(0,1.5e5,5e4)
axis(side=1,at=xticks,label=as.expression(sciNotation(xticks, 1)))
axis(side=2,at=axTicks(2),las=2,label=as.expression(sciNotation(axTicks(2), 1)))
box()

a_m0_fit = lm(y ~ 0 + x)
abline(a_m0_fit,lty=2)
Rsq = as.character(round(summary(a_m0_fit)$r.squared,2))
legend("topleft",as.expression(bquote(R^2 == .(Rsq))),
       cex=0.8,bty='n')
if (save_eps) dev.off()