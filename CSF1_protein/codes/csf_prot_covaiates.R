## csf1_prot_covariates.R
## covariates of csf1  proteomics
if(T){
  rm(list = ls())
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF1/")
  
  ## load data ####
  load("./codes/CSF1_proteomics_normalized.Rdata")
  
  ## plot covariates ####
  cols = c("grey","green")
  sex = unique(datMeta$Sex)
  for(s in sex){
    
    pdf(paste0("./results/figures/covariates/csf1_prot_covaritates_",s,".pdf"), 
        height = 5, 
        width =6,
        title = paste0(s, "covariates"))
    par(mfrow = c(2,3))
    
    #---1. Subjects
    plot(datMeta$Diagnosis, col=cols, main= "conditions")
    
    # #---2. Sex
    # A = chisq.test(as.factor(datMeta$Sex), datMeta$Diagnosis)
    # p = A$p.value
    # plot(as.factor(datMeta$Sex) ~ datMeta$Diagnosis, col= c("blue", "magenta"), 
    #      main=paste("Sex \np=", signif(p,2)), ylab="", xlab="")
    # 
    #---3. Age
    A= anova(lm(as.numeric(datMeta$Age) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$Age ~ datMeta$Diagnosis, col= cols, 
         main = paste("Age \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    #---4. MoCA
    A= anova(lm(as.numeric(datMeta$MoCA) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$MoCA ~ datMeta$Diagnosis, col= cols, 
         main = paste("MoCA \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    #---5. pTau
    A= anova(lm(as.numeric(datMeta$pTau.ELISA) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$pTau.ELISA ~ datMeta$Diagnosis, col= cols, 
         main = paste("pTau \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    #---6. AB42
    A= anova(lm(as.numeric(datMeta$AB42.ELISA) ~ datMeta$Diagnosis))
    p= A$"Pr(>F)"[1]
    plot(datMeta$AB42.ELISA ~ datMeta$Diagnosis, col= cols, 
         main = paste("AB42 \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    
    #---7 APOE
    A = chisq.test(as.factor(datMeta$APOE.Genotype), datMeta$Diagnosis)
    p = A$p.value
    plot(as.factor(datMeta$APOE.Genotype) ~ datMeta$Diagnosis, 
         col= cols, 
         main = paste("APOE genotype \np=", signif(p,2), sep=""), ylab="", xlab="")
    
    
    dev.off()
  }
}