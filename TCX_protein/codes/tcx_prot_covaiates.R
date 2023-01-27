## tcx_prot_covariates.R
## covariates of temporal cortex proteomics
if(T){
  rm(list = ls())
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_temoral/")

## load data ####
load("./codes/tcx_proteomics_normalized.Rdata")

## plot covariates ####
cols = c("grey","blue")
sex = unique(datMeta$Sex)

for(s in sex){
  id = which(datMeta$Sex==s)
  meta = datMeta[id,]
  
  pdf(paste("./results/figures/covariates/tcx_prot_covaritates_",s,".pdf",sep=""), 
      height = 6, 
      width =6)
  par(mfrow = c(2,4))
  
  #---1. Subjects
  plot(meta$Diagnosis, col=cols, main= "conditions")
  
  #---2. Age
  A= anova(lm(as.numeric(meta$AgeAtDeath) ~ meta$Diagnosis))
  p= A$"Pr(>F)"[1]
  plot(meta$AgeAtDeath ~ meta$Diagnosis, col= cols, 
       main = paste(s,"- Age \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  ##--4 PMI
  A= anova(lm(as.numeric(meta$PMI) ~ meta$Diagnosis))
  p= A$"Pr(>F)"[1]
  plot(meta$PMI ~ meta$Diagnosis, 
       col= cols, 
       main = paste(s,"- PMI \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  #---5. RIN
  A= anova(lm(as.numeric(meta$RIN) ~ meta$Diagnosis))
  p= A$"Pr(>F)"[1]
  plot(meta$RIN ~ meta$Diagnosis, 
       col= cols, 
       main = paste(s,"- RIN \np=", signif(p,2), sep=""), ylab="", xlab="")
  #---5. Braak
  A= anova(lm(as.numeric(meta$Braak) ~ meta$Diagnosis))
  p= A$"Pr(>F)"[1]
  plot(meta$Braak ~ meta$Diagnosis, 
       col= cols, 
       main = paste(s,"- Braak \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  #---6 APOE
  A = chisq.test(as.factor(meta$ApoE.Genotype), meta$Diagnosis)
  p = A$p.value
  plot(as.factor(meta$ApoE.Genotype) ~ meta$Diagnosis, 
       col= cols, 
       main = paste(s,"- APOE genotype \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  
  dev.off()
}
}