## csf2_prot_covaiates.R
## covariates of csf2
if(T){ 
 rm(list = ls())
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF2/")
  
  ## load data ####
  load("./codes/CSF2_proteomics_normalized.Rdata")
  datMeta$ptau = as.numeric(datMeta$ptau)
  datMeta$ptau[is.na(datMeta$ptau)]= mean(datMeta$ptau,na.rm=T)
  datMeta$Moca[is.na(datMeta$Moca)]= mean(datMeta$Moca,na.rm=T)
  datMeta$AB.TAUratio[is.na(datMeta$AB.TAUratio)]= mean(datMeta$AB.TAUratio,na.rm=T)
  
  ## plot covariates ####
  cols = c("grey","blue")
  sex = unique(datMeta$Sex)
  for(s in sex){
    
    pdf(paste0("./results/figures/covariates/csf2_prot_covaritates_",s,".pdf"), 
        height = 5, 
        width =6)
    par(mfrow = c(2,3))
    
    #---1. Subjects
    plot(datMeta$Diagnosis, col=cols, main= paste0(s,"- conditions"))
    
    #-- Age
    A= anova(lm(as.numeric(datMeta$Age) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$Age ~ datMeta$Diagnosis, col= cols, 
         main = paste(s,"-Age \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    #-- Moca
    A= anova(lm(as.numeric(datMeta$Moca) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$Moca ~ datMeta$Diagnosis, col= cols, 
         main = paste(s,"-Moca \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    #- pTau
    A= anova(lm(as.numeric(datMeta$ptau) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$ptau ~ datMeta$Diagnosis, col= cols, 
         main = paste(s,"-pTau \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    # AB42
    A= anova(lm(as.numeric(datMeta$Ab42) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$Ab42 ~ datMeta$Diagnosis, col= cols, 
         main = paste(s,"-AB42 \np=", signif(p,2), sep=""), ylab="", xlab="")
    # AB/tau ratio
    A= anova(lm(as.numeric(datMeta$AB.TAUratio) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$AB.TAUratio ~ datMeta$Diagnosis, col= cols, 
         main = paste(s,"-AB42/tau ratio \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    
    
    dev.off()
  }
}