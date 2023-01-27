## csf_prot_covaiates.R
if(T){
rm(list=ls())
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF_bader/")

## load expression data
datExp = read.csv("../../data/proteomics/bader_et_al/CSF/datexp.csv")

gene =str_split(datExp$PG.Genes,pattern = "-",simplify = T)[,1]
id = !duplicated(gene)
gene = gene[id]
datExp = datExp[id,]
rownames(datExp) = gene
datExp = datExp[,-c(1:3)]
colnames(datExp) = gsub("X*\\..","",colnames(datExp))

## covariates
datMeta = read.csv("../../data/proteomics/bader_et_al/CSF/metadata.csv")
datMeta = as.data.frame(datMeta)
id = which(datMeta$collection.site == "Sweden")
datMeta = datMeta[id,]
datExp = datExp[,id]

datMeta$sample.name = gsub("sample.*","sample",datMeta$sample.name)
colnames(datExp) = gsub("sample.*","sample",colnames(datExp))

id = match(datMeta$sample.name, colnames(datExp))
## keep only AD and CTL
datMeta$Dx = gsub("biochemical ","",datMeta$biochemical.AD.classification)
datMeta$Dx = gsub("control","CTL",datMeta$Dx)
id = which(datMeta$Dx %in% c("CTL","AD"))
datMeta = datMeta[id,]
datExp = datExp[,id]
datMeta$Dx = factor(datMeta$Dx, levels = c("CTL","AD"))

## match with expression data 
rownames(datMeta) = datMeta$Sample
id = match(rownames(datMeta),colnames(datExp))
datExp = datExp[,id]
all(rownames(datMeta)==colnames(datExp))

datMeta$Diagnosis = factor(datMeta$Group, levels = c("CTL","AD"))

## impute missing values ####
# exp = impute::impute.knn(as.matrix(datExp),k = 5, rowmax = 0.5,
#                  colmax = 0.5, maxp = 500, rng.seed=12)
exp = imputeLCMD::impute.QRILC(as.matrix(datExp))
datExp = exp[[1]]

save(file = "./codes/CSF1_proteomics_normalized.Rdata",datExp,datMeta,attr)



## load data ####
load("./codes/CSF1_proteomics_normalized.Rdata")
colnames(datMeta)

## plot covariates ####
condition = T
cols = c("grey","green")

if(condition){
  
  pdf("./results/figures/covariates/csf_prot_covaritates.pdf", 
      height = 5, 
      width =6)
  par(mfrow = c(2,4))
  
  #---1. Subjects
  plot(datMeta$Diagnosis, col=cols, main= "conditions")
  
  #---2. Sex
  A = chisq.test(as.factor(datMeta$Sex), datMeta$Diagnosis)
  p = A$p.value
  plot(as.factor(datMeta$Sex) ~ datMeta$Diagnosis, col= c("blue", "magenta"), 
       main=paste("Sex \np=", signif(p,2)), ylab="", xlab="")
  
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
