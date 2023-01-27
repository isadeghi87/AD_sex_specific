###--msbb_covariates.R
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/bulk_analysis/replication_msbb/")

#libraries

## load data ####
load("./codes/MSBB_normalized.Rdata")
colnames(datMeta)
## plot covariates ####
condition = T
cols = c("grey","orange")

if(condition){
  
  pdf("./results/figures/covariates/msbb_covaritates.pdf", 
      height = 5, 
      width =8)
  par(mfrow = c(2,5))
  
  #---1. Subjects
  plot(datMeta$Dx, col=cols, main= "conditions")
  
  #---2. Sex
  A = chisq.test(datMeta$Sex, datMeta$Dx)
  p = A$p.value
  plot(datMeta$Sex ~ datMeta$Dx, col= c("blue", "magenta"), 
       main=paste("Sex \np=", signif(p,2)), ylab="", xlab="")
  
  #---3. Age
  A= anova(lm(as.numeric(datMeta$Age) ~ datMeta$Dx))
  p= A$"Pr(>F)"[1]
  plot(datMeta$Age ~ datMeta$Dx, col= cols, 
       main = paste("Age \np=", signif(p,2), sep=""), ylab="year", xlab="")
  
  ##--4 PMI
  A= anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Dx))
  p= A$"Pr(>F)"[1]
  plot(datMeta$PMI ~ datMeta$Dx, 
       col= cols, 
       main = paste("PMI \np=", signif(p,2), sep=""), ylab="min", xlab="")
  
  #---5. RIN
  A= anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Dx))
  p= A$"Pr(>F)"[1]
  plot(datMeta$RIN ~ datMeta$Dx, 
       col= cols, 
       main = paste("RIN \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  #---6. CDR
  A= anova(lm(as.numeric(datMeta$CDR) ~ datMeta$Dx))
  p= A$"Pr(>F)"[1]
  plot(datMeta$CDR ~ datMeta$Dx, 
       col= cols, 
       main = paste("CDR \np=", signif(p,2), sep=""), ylab="", xlab="")
  #---7. plaque mean
  A= anova(lm(as.numeric(datMeta$PlaqueMean) ~ datMeta$Dx))
  p= A$"Pr(>F)"[1]
  plot(datMeta$PlaqueMean ~ datMeta$Dx, 
       col= cols, 
       main = paste("Plaque mean \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  #---8. NP
  A= anova(lm(as.numeric(datMeta$NP.1) ~ datMeta$Dx))
  p= A$"Pr(>F)"[1]
  plot(datMeta$NP.1 ~ datMeta$Dx, 
       col= cols, 
       main = paste("NP \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  #---9. number of reads
  A = anova(lm(as.numeric(datMeta$TotalReads) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   
  plot(as.numeric(datMeta$TotalReads) ~ datMeta$Dx, col=cols,
       main=paste("Reads, p=", signif(p,2)), ylab="", xlab="")
  #---10. number_of_mapped_reads
  A = anova(lm(as.numeric(datMeta$Mapped) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   
  plot(as.numeric(datMeta$Mapped) ~ datMeta$Dx,col=cols,
       main=paste("Mapped Reads, p=", signif(p,2)), ylab="", xlab="")

  
  dev.off()
}
