### csf_proteomics_differential.R
## we compute differential expression for CSF2 proteomics data 
## for AsymAD obtained from lenora et al
if(T){
  rm(list=ls())
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF2/")

library(reshape); library(readxl)
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);library(tidyr)
library(statmod);library(edgeR);library(MSnbase)

## load expression data
datExp = read.csv("../../data/proteomics/Erik_et_al/CSF2/CSF2_unregressed_Batch-corrected_cleanDat.csv",
                  row.names = 1)

gene =str_split(rownames(datExp),pattern = "\\|",simplify = T)[,1]
# protein_ID =str_split(rownames(datExp),pattern = "\\|",simplify = T)[,2]
# 
# 
# id = which(duplicated(gene))
# gene = gene[-id]
# datExp = datExp[-id,]
# rownames(datExp) = gene

## covariates
datMeta = read.csv("../../data/proteomics/Erik_et_al/CSF2/CSF2_Traits.csv")

## keep only AD and control
datMeta$ClinicalGroup = gsub("NORMAL","CTL",datMeta$ClinicalGroup)

## keep only AD and CTL
datMeta = subset(datMeta,ClinicalGroup %in% c("CTL","MCI"))

datMeta$Sex[datMeta$Sex==1] = "male"
datMeta$Sex[datMeta$Sex==0] = "female"

## match with expression data 
rownames(datMeta) = datMeta$batch.channel
id = match(rownames(datMeta),colnames(datExp))
datExp = datExp[,id]
all(rownames(datMeta)==colnames(datExp))

datMeta$Diagnosis = factor(datMeta$ClinicalGroup, levels = c("CTL","MCI"))

## impute missing values ####
# exp = impute::impute.knn(as.matrix(datExp),k = 5, rowmax = 0.5,
#                  colmax = 0.5, maxp = 500, rng.seed=12)
exp = imputeLCMD::impute.QRILC(as.matrix(datExp))
datExp = exp[[1]]

save(file = "./codes/CSF2_proteomics_normalized.Rdata",datExp,datMeta,gene)

sex = unique(datMeta$Sex)
## differential expression
stats = data.frame()
for(s in sex){
  id = which(datMeta$Sex == s)
  meta = datMeta[id,]
  exp = datExp[,id]
  mod = model.matrix( ~ Diagnosis,data = meta)
  
  fit <- lmFit(exp, mod)
  efit <- eBayes(fit, trend = T, robust = T)
  meta <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
  meta <- meta[,c(1,4,5)]
  colnames(meta) = c("logFC", "p.value", "fdr")
  sig = subset(meta, p.value<0.05) 
  sig$sex = s
  stats = rbind(stats,sig)
  
}

write.table(stats,"./results/tables/differential/CSF2_proteomics_sex_stats.txt",
            col.names = T,row.names = T)
}