# AD_PA_sex_covariates.R
#--- meta analysis for AD and PA -sex-specific
## here we analyze only temporal cortex
## data is the same as for bulk analysis
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/")
load("./AD_PA_CTL/bulk_analysis/codes/covariates/ad_pa_general.Rdata")


## plot covariates for each sex ####
sex = unique(datmeta$Sex)

for(s in sex){
  idx = which(datmeta$Sex == s)
  meta = datmeta[idx,]
  exp = datExp[,idx]
  
  pdf(paste("./AD_PA_CTL/sex_specific/results/figures/covariates/covaritates_",s,".pdf",sep=""), 
      height = 7, 
      width =10)
  par(mfrow = c(3,4))
  
  #---1. Subjects
  plot(meta$Dx, col=c("grey", "yellow","black"), main= "conditions")
  
  #---2. Age
  A= anova(lm(as.numeric(meta$Age) ~ meta$Dx))
  p= A$"Pr(>F)"[1]
  plot(meta$Age ~ meta$Dx, col= c("grey", "yellow","black"), 
       main = paste("Age \np=", signif(p,2), sep=""), ylab="year", xlab="")
  
  #---3. APOE
  A= anova(lm(as.numeric(meta$APOE) ~ meta$Dx))
  p= A$"Pr(>F)"[1]
  plot(meta$Dx ~ meta$APOE , col= c("black", "yellow", "grey"), 
       main = paste("APOE \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  #---4. RIN
  A= anova(lm(as.numeric(meta$RIN) ~ meta$Dx))
  p= A$"Pr(>F)"[1]
  plot(meta$RIN ~ meta$Dx, col= c("grey", "yellow","black"), 
       main = paste("RIN \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  
  #---5. number of reads
  A = anova(lm(as.numeric(meta$Number_of_reads) ~ meta$Dx)); p = A$"Pr(>F)"[1];   
  plot(as.numeric(meta$Number_of_reads) ~ meta$Dx, col=c("grey", "yellow", "black"),
       main=paste("Reads, p=", signif(p,2)), ylab="", xlab="")
  #---6. number_of_mapped_reads
  A = anova(lm(as.numeric(meta$Numer_of_mapped_reads) ~ meta$Dx)); p = A$"Pr(>F)"[1];   
  plot(as.numeric(meta$Numer_of_mapped_reads) ~ meta$Dx,col=c("grey", "yellow", "black"),
       main=paste("Mapped Reads, p=", signif(p,2)), ylab="", xlab="")
  #----7. GC_percentage
  A = anova(lm(as.numeric(meta$GC_percentage) ~ meta$Dx)); p = A$"Pr(>F)"[1];   
  plot(as.numeric(meta$GC_percentage) ~ meta$Dx,col=c("grey", "yellow", "black"),
       main=paste("GC_Percentage, p=", signif(p,2)), ylab="", xlab="")
  #---8-11. SeqPC
  for(pc in paste("seqPC",1:4,sep=""))  {
    A = anova(lm(as.numeric((meta[,pc])) ~ meta$Dx));   
    p = A$"Pr(>F)"[1];   
    plot((meta[,pc]) ~ meta$Dx, col=c("grey", "yellow", "black"),
         main=paste(pc, ", p=", signif(p,2)), ylab="", xlab="")
  }
  
  
  dev.off()
}
