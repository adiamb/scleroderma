library(ineq)
library(tidyr)
library(alakazam)
library(data.table)
chao1 = function(abunds) {
  f1 = sum(abunds==1)
  f2 = sum(abunds==2)
  D = length(abunds) + f1^2/(2*f2)
  return(D)
}

save_eps=T
calculate=T

armA_pats = ####
armB_pats = ####
healthy_pats = ####
id_to_dummy =  c(paste("SR",1:length(armA_pats),sep=""),
                 paste("SP",1:length(armB_pats),sep=""),
                 paste("H",1:length(healthy_pats),sep=""))
names(id_to_dummy) = c(armA_pats,armB_pats,healthy_pats)
dummy_to_id = c(armA_pats,armB_pats,healthy_pats)
names(dummy_to_id) = c(paste("SR",1:length(armA_pats),sep=""),
                       paste("SP",1:length(armB_pats),sep=""),
                       paste("H",1:length(healthy_pats),sep=""))

pats = ####
path = ####
fig_path = ####
mut_field = "NMUT_GERMLINE" 

if (calculate) {
s0bl = m0bl = m1bl = m0 = dm0dt = t0 = m0pre = t0pre = list()
for (pat in pats) {
  file = paste(path,pat,"_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab",sep="")
  data = readChangeoDb(file)
  data$ISOTYPE = gsub("1|2|3|4","",data$PRIMER)

  M0bl = chao1(subset(data,WEEK==0 & ISOTYPE=="IgM" & eval(parse(text=paste(mut_field,"==0",sep=""))),
                    DUPCOUNT,drop=T))
  M1bl = chao1(subset(data,WEEK==0 & ISOTYPE=="IgM" & eval(parse(text=paste(mut_field,"==1",sep=""))),
                      DUPCOUNT,drop=T))
  S0bl = chao1(subset(data,WEEK==0 & !(ISOTYPE %in% c("IgD","IgM")) & 
                      eval(parse(text=paste(mut_field,"==0",sep=""))),DUPCOUNT,drop=T))

  M0 = aggregate(DUPCOUNT ~ WEEK_CORRECTED+WEEK,FUN=chao1,
                 data=subset(data,ISOTYPE=="IgM" & eval(parse(text=paste(mut_field,"==0",sep="")))))
  M0 = M0[order(M0$WEEK_CORRECTED),]
  
  M0pre = subset(M0, WEEK < 24)
  T0pre = M0pre$WEEK_CORRECTED
  M0pre = M0pre$DUPCOUNT
  
  M0 = subset(M0, WEEK >= 24)
  
  #*#*#* smoothing (window of size 2)
  T0 = M0$WEEK_CORRECTED[1:(nrow(M0)-1)]/2+M0$WEEK_CORRECTED[2:nrow(M0)]/2
  M0 = M0$DUPCOUNT[1:(nrow(M0)-1)]/2+M0$DUPCOUNT[2:nrow(M0)]/2
  #*#*#*
  
  #*#*#* estimate derivative
  dM0dt = diff(M0)/diff(T0)
  M0 = (M0[1:(length(M0)-1)]+M0[2:length(M0)])/2
  T0 = (T0[1:(length(T0)-1)]+T0[2:length(T0)])/2
  #*#*#*
  
  #*#*#* remove NAs
  any_NAs = (is.na(M0) | is.na(dM0dt))
  M0 = M0[!any_NAs]
  dM0dt = dM0dt[!any_NAs]
  T0 = T0[!any_NAs]
  #*#*#*
  

  plot(M0,dM0dt,main=pat)
  m0[[pat]] = M0
  m0pre[[pat]] = M0pre
  t0pre[[pat]] = T0pre
  dm0dt[[pat]] = dM0dt
  t0[[pat]] = T0
  s0bl[[pat]] = S0bl
  m0bl[[pat]] = M0bl
  m1bl[[pat]] = M1bl
}
saveRDS(m0,paste(path,"m0","_chao_midpointDerivative.RDS",sep=""))
saveRDS(m0pre,paste(path,"m0pre","_chao_midpointDerivative.RDS",sep=""))
saveRDS(t0pre,paste(path,"t0pre","_chao_midpointDerivative.RDS",sep=""))
saveRDS(dm0dt,paste(path,"dm0dt","_chao_midpointDerivative.RDS",sep=""))
saveRDS(t0,paste(path,"t0","_chao_midpointDerivative.RDS",sep=""))
saveRDS(s0bl,paste(path,"s0bl","_chao_midpointDerivative.RDS",sep=""))
saveRDS(m1bl,paste(path,"m1bl","_chao_midpointDerivative.RDS",sep=""))
saveRDS(m0bl,paste(path,"m0bl","_chao_midpointDerivative.RDS",sep=""))
}


m0 = readRDS(paste(path,"m0","_chao_midpointDerivative.RDS",sep=""))
m0pre = readRDS(paste(path,"m0pre","_chao_midpointDerivative.RDS",sep=""))
t0pre = readRDS(paste(path,"t0pre","_chao_midpointDerivative.RDS",sep=""))
dm0dt = readRDS(paste(path,"dm0dt","_chao_midpointDerivative.RDS",sep=""))
t0 = readRDS(paste(path,"t0","_chao_midpointDerivative.RDS",sep=""))
s0bl = readRDS(paste(path,"s0bl","_chao_midpointDerivative.RDS",sep=""))
m1bl = readRDS(paste(path,"m1bl","_chao_midpointDerivative.RDS",sep=""))
m0bl = readRDS(paste(path,"m0bl","_chao_midpointDerivative.RDS",sep=""))

### Fit linear models jointly and test if coefficients are the same (http://stats.stackexchange.com/questions/12797/how-to-test-whether-a-regression-coefficient-is-moderated-by-a-grouping-variable)
joint_df = data.frame()
for (pat in names(m0)) {
  M0 = m0[[pat]]
  dM0dt = dm0dt[[pat]]
  T0 = t0[[pat]]
  #*#*#* remove earliest point for SR1 (only interested in monotonically increasing tail of T0-vs-M0 curve)
  if (pat==###) {
    cond = (T0 >= 37)
    M0 = M0[cond]
    dM0dt = dM0dt[cond]
  }
  #*#*#*
  joint_df = rbind(joint_df,data.frame(x=M0,y=dM0dt,pat=pat))
}
joint_df$x2 = joint_df$x
joint_df$pat2 = paste("x",joint_df$pat,sep="")
joint_df = spread(joint_df,key=pat2,value=x2,fill=0)
joint_df = spread(cbind(joint_df,dummy=1),key=pat,value=dummy,fill=0)
reference_pat = ###SR6
joint_df = joint_df[,!grepl(reference_pat,colnames(joint_df))]
joint_fit = lm(y ~ ., joint_df)
summary(joint_fit)
# save parameter results
a_jointfit = delta_jointfit = numeric()
for (pat in pats) {
  if (pat==reference_pat) {
    a_jointfit[pat] = coef(joint_fit)["(Intercept)"]
    delta_jointfit[pat] = -coef(joint_fit)["x"]
  } else {
    a_jointfit[pat] = coef(joint_fit)["(Intercept)"] + coef(joint_fit)[pat]
    delta_jointfit[pat] = -coef(joint_fit)["x"] - coef(joint_fit)[paste("x",pat,sep="")]
  }
}
a_jointfit = na.omit(a_jointfit)
saveRDS(a_jointfit,paste(path,"a_jointfreefit","_chao.RDS",sep=""))
saveRDS(delta_jointfit,paste(path,"delta_jointfreefit","_chao.RDS",sep=""))
### Make table
df = as.data.frame(summary(joint_fit)$coefficients[,c("Estimate","Std. Error","Pr(>|t|)")])
setnames(df,"Pr(>|t|)","p-value")
rownames(df)[rownames(df)=="(Intercept)"] = paste("intercept(",id_to_dummy[reference_pat],")",sep="")
rownames(df)[rownames(df)=="x"] = paste("slope(",id_to_dummy[reference_pat],")",sep="")
pats=####
for (pat in setdiff(pats,reference_pat)) {
  rownames(df)[rownames(df)==pat] = paste("intercept(",id_to_dummy[pat],") - ",
                                          "intercept(",id_to_dummy[reference_pat],")",sep="")
  rownames(df)[rownames(df)==paste("x",pat,sep="")] = paste("slope(",id_to_dummy[pat],") - ",
                                          "slope(",id_to_dummy[reference_pat],")",sep="")
}
df = df[order(rownames(df)),]
View(df)

### Fit linear model, requiring slope to be same for all patients,
### plot and save parameters
joint_df = data.frame()
for (pat in names(m0)) {
  M0 = m0[[pat]]
  dM0dt = dm0dt[[pat]]
  T0 = t0[[pat]]
  #*#*#* remove earliest point for SR1 (only interested in monotonically increasing tail of T0-vs-M0 curve)
  if (pat==###) {
    cond = (T0 >= 37)
    M0 = M0[cond]
    dM0dt = dM0dt[cond]
  }
  #*#*#*
  joint_df = rbind(joint_df,data.frame(x=M0,y=dM0dt,pat=pat))
}
joint_df = spread(cbind(joint_df,dummy=1),key=pat,value=dummy,fill=0)
reference_pat = ###SR6
joint_df = joint_df[,!grepl(reference_pat,colnames(joint_df))]
joint_fit = lm(y ~ ., joint_df)
summary(joint_fit)
# save resulting parameters
delta_jointfit = -coef(joint_fit)["x"]
a_jointfit = b_jointfit = c1_jointfit = c2_jointfit = numeric()
a_error = b_error = c1_error = c2_error = numeric()
jointfit_errors = coef(summary(joint_fit))[,"Std. Error"]
delta_error = jointfit_errors["x"]
for (pat in pats) {
  if (pat==reference_pat) {
    a_jointfit[pat] = coef(joint_fit)["(Intercept)"]
    a_error[pat] = jointfit_errors["(Intercept)"]
  } else {
    a_jointfit[pat] = coef(joint_fit)["(Intercept)"] + coef(joint_fit)[pat]
    a_error[pat] = sqrt((jointfit_errors["(Intercept)"])^2 + (jointfit_errors[pat])^2)
  }
  S0bl = s0bl[[pat]] 
  M0bl = m0bl[[pat]] 
  M1bl = m1bl[[pat]]
  b_jointfit[pat] = delta_jointfit*M1bl/M0bl
  c1_jointfit[pat] = delta_jointfit*S0bl/(M0bl+S0bl)
  c2_jointfit[pat] = delta_jointfit-b_jointfit[pat]-c1_jointfit[pat] 
  b_error[pat] = delta_error*M1bl/M0bl
  c1_error[pat] = delta_error*S0bl/(M0bl+S0bl)
  c2_error[pat] = sqrt(delta_error^2+(b_error[pat])^2+(c1_error[pat])^2) 
}
a_jointfit = na.omit(a_jointfit)
saveRDS(a_jointfit,paste(path,"a_constrainedfit","_chao.RDS",sep=""))
saveRDS(delta_jointfit,paste(path,"delta_constrainedfit","_chao.RDS",sep=""))
saveRDS(b_jointfit,paste(path,"b_constrainedfit","_chao.RDS",sep=""))
saveRDS(c1_jointfit,paste(path,"c1_constrainedfit","_chao.RDS",sep=""))
saveRDS(c2_jointfit,paste(path,"c2_constrainedfit","_chao.RDS",sep=""))
saveRDS(a_error,paste(path,"a_error_constrainedfit","_chao.RDS",sep=""))
saveRDS(delta_error,paste(path,"delta_error_constrainedfit","_chao.RDS",sep=""))
saveRDS(b_error,paste(path,"b_error_constrainedfit","_chao.RDS",sep=""))
saveRDS(c1_error,paste(path,"c1_error_constrainedfit","_chao.RDS",sep=""))
saveRDS(c2_error,paste(path,"c2_error_constrainedfit","_chao.RDS",sep=""))

# Fig 4B
df = data.frame()
for (pat in names(a_jointfit)) {
  df[id_to_dummy[pat],"a"] = paste(a_jointfit[pat],a_error[pat],sep=" +/-")
  df[id_to_dummy[pat],"b"] = paste(b_jointfit[pat],b_error[pat],sep=" +/-")
  df[id_to_dummy[pat],"c1"] = paste(c1_jointfit[pat],c1_error[pat],sep=" +/-")
  df[id_to_dummy[pat],"c2"] = paste(c2_jointfit[pat],c2_error[pat],sep=" +/-")
}
View(df)

# Fig 4A
delta = delta_jointfit
fig4atable = data.frame()
for (pat in names(a_jointfit)) {
  if (save_eps) {
    pdf(file=paste(fig_path,pat,"-fit-exponential.eps",sep=""), 
        width=3,height=2.5,paper="a4",
        colormodel="rgb",pointsize = 10)
  }
  par(mar=c(6,6,2,2))
  
  A = a_jointfit[pat]
  M0pre = m0pre[[pat]]
  T0pre = t0pre[[pat]]
  M0 = m0[[pat]]
  T0 = t0[[pat]]
  #*#*#* remove earliest point for SR1 (only interested in monotonically increasing tail of T0-vs-M0 curve)
  if (pat==###) {
    cond = (T0 >= 37)
    M0 = M0[cond]
    T0 = T0[cond]
  }
  #*#*#*
  T0_pred = seq(T0[1],T0[length(T0)],1)
  M0_pred = A/delta - (A/delta - M0[1])*exp(-delta*(T0_pred-T0_pred[1]))
  plot(T0,M0,ylim=range(c(M0,M0_pred)),
       xlab="Week",ylab=expression(M[0]),
       ann=F,axes=F,pch=20)
  fig4atable = rbind(fig4atable,data.frame(Participant=rep(id_to_dummy[pat],length(T0)),
                                           Week=T0,M0=M0))
  lines(T0_pred,M0_pred,lty=2)
  axis(side=1,at=axTicks(1))
  mtext(side=1,line=2,"Week")
  axis(side=2,at=axTicks(2),las=2,label=as.expression(sciNotation(axTicks(2), 1)))
  mtext(side=2,line=4,expression(M[0]))
  box()
  legend("bottomright",legend=id_to_dummy[pat],bty='n')
  if (save_eps) dev.off()
}