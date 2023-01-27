### temporal_proteomics_differential.R
## we compute differential expression for temporal proteomics data 
## for Alzheimer obtained from Eric et al, Nat Med
if(T){
  rm(list = ls())
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_temoral/")

library(reshape); 
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);library(tidyr)
library(statmod);library(edgeR)

## load expression data
exp = read.csv("../../data/proteomics/Erik_et_al/temporal/TCX_proteomics_data.csv",
               row.names = 1)

## covariates
datMeta = read.csv("../../data/proteomics/Erik_et_al/temporal/Traits.csv")

## keep only AD and control
datMeta$Diagnosis = gsub("Control","CTL",datMeta$Diagnosis)
datMeta = subset(datMeta,Diagnosis %in% c("AD","CTL"))
rownames(datMeta) = datMeta$SampleID
datMeta$Diagnosis = factor(datMeta$Diagnosis, levels = c("CTL","AD"))

## match with expression data 
id = match(colnames(exp),rownames(datMeta))
datMeta = datMeta[id,]
all(rownames(datMeta)==colnames(exp))

## log transformation
design = model.matrix(~Diagnosis,datMeta)

## handle NAs ## we try to use mean of counts
for(i in 1:nrow(exp)){
  exp[i,][is.na(exp[i,])]= mean(as.numeric(exp[i,]),na.rm = T)
}

v = voom(exp,design)
datExp = v$E
datMeta$Sex[datMeta$Sex==1] = "male"
datMeta$Sex[datMeta$Sex==0] = "female"
save(file = "./codes/tcx_proteomics_normalized.Rdata",datExp,datMeta)

sumstats = data.frame()

for (s in unique(datMeta$Sex)){
## differential expression
mod = model.matrix(~ Diagnosis + AgeAtDeath + PMI ,data = datMeta)

fit <- lmFit(datExp, mod)
efit <- eBayes(fit, trend = T, robust = T)
meta <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
meta <- meta[,c(1,4,5)]
colnames(meta) = c("logFC", "p.value", "fdr")
meta$gene = rownames(meta)
rownames(meta) = NULL
meta$sex = s
sumstats = rbind(sumstats,meta)
}
write.table(sumstats,"./results/tables/differential/temporal_proteomics_sex_sumstats.txt",col.names = T,row.names = T)
}
