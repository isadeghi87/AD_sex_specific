### temporal_proteomics_differential.R
## we compute differential expression for temporal proteomics data 
## for Alzheimer obtained from Eric et al, Nat Med
if(T){
  rm(list = ls())
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")

library(reshape); 
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);library(tidyr)
library(statmod);library(edgeR)

## load expression data
load("./data/normalized.Rdata")


## handle NAs ## we try to use mean of counts
for(i in 1:nrow(exp)){
  exp[i,][is.na(exp[i,])]= mean(as.numeric(exp[i,]),na.rm = T)
}

v = voom(exp,design)
datExp = v$E
datMeta$Sex[datMeta$Sex==1] = "male"
datMeta$Sex[datMeta$Sex==0] = "female"

sumstats = data.frame()
for (s in unique(datMeta$Sex)){
## differential expression
mod = model.matrix(~ Diagnosis,data = datMeta)


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
