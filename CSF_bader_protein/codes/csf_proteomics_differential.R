### csf_proteomics_differential.R
## we compute differential expression for CSF1 proteomics data 
## for AD and AsymAD obtained from lenora et al

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF_bader//")

library(reshape); library(readxl)
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);library(tidyr)
library(statmod);library(edgeR);library(MSnbase)

## load expression data
datExp = read.csv("../../data/proteomics/Erik_et_al/CSF1/2b.Unregressed_Batch-corrected_cleanDat.csv",
                  row.names = 1)

gene =str_split(rownames(datExp),pattern = "\\|",simplify = T)[,1]
protein_ID =str_split(rownames(datExp),pattern = "\\|",simplify = T)[,2]

### annotate 
getinfo <- c( "ensembl_gene_id",
              "external_gene_name")

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")
attr = getBM(
  attributes = getinfo,
  filters = "external_gene_name",
  values = unique(gene),
  mart = mart)

id = duplicated(attr$external_gene_name)
attr = attr[!id,]
id = match(attr$external_gene_name,gene)
datExp = datExp[id,]
rownames(datExp) = attr$ensembl_gene_id

## covariates
datMeta = read.csv("../../data/proteomics/Erik_et_al/CSF1/CSF1_Traits.csv")
datMeta = as.data.frame(datMeta)

## keep only AD and control
datMeta$Group = gsub("Control","CTL",datMeta$Group)

## keep only AD and CTL
datMeta = subset(datMeta,Group %in% c("CTL","AD"))

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

## differential expression
mod = model.matrix(~ Diagnosis  ,data = datMeta)

fit <- lmFit(datExp, mod)
efit <- eBayes(fit, trend = T, robust = T)
meta <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
meta <- meta[,c(1,3,4,5)]
colnames(meta) = c("logFC", "t", "p.value", "fdr")
sig_proteins = subset(meta, p.value<0.05) 

sig_proteins$gene =attr$external_gene_name[match(rownames(sig_proteins),attr$ensembl_gene_id)]

write.table(sig_proteins,"./results/tables/differential/CSF1_proteomics_DE.txt",
            col.names = T,row.names = T)
