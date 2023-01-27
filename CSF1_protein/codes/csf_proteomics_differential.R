### csf_proteomics_differential.R
## we compute differential expression for CSF1 proteomics data 
## for AD and AsymAD obtained from lenora et al
if(T){
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF1/")

library(reshape); library(readxl);library(WGCNA)
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);library(tidyr)
library(statmod);library(edgeR);library(MSnbase)

## load expression data
datExp = read.csv("../../data/proteomics/Erik_et_al/CSF1/2b.Unregressed_Batch-corrected_cleanDat.csv",
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
datMeta = read.csv("../../data/proteomics/Erik_et_al/CSF1/CSF1_Traits.csv")
datMeta = as.data.frame(datMeta)
datMeta$Sex[datMeta$Sex==1] = "male"
datMeta$Sex[datMeta$Sex==0] = "female"
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

## check PCA for outliers
pc = prcomp(t(datExp))

p = factoextra::fviz_pca_ind(pc,
                             # label="none",
                             col.ind  = meta$Diagnosis)+
  labs(title = s)+
  scale_color_brewer(palette="Set1")+ 
  theme_bw()+
  theme(panel.grid = element_blank())

## b37.131N is an outlier
id = which(datMeta$SampleID!="b37.131N")
datMeta = datMeta[id,]
datExp = datExp[,id]

## check more outliers
normadj <- (0.5+0.5*bicor(datExp))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,
     col=c("blue" , "red"),
     pch=19, xlab = "",ylab= "Z-score")
legend("bottomleft",legend = unique(datMeta$Diagnosis), 
       col = c("blue", "red"),pch=19,cex=.7)
abline(h=-2, lty=2)
keep = (z.ku >= -2)
datMeta = datMeta[keep,]
datExp = datExp[,keep]

save(file = "./codes/CSF1_proteomics_normalized.Rdata",datExp,datMeta,gene)

sex = unique(datMeta$Sex)
## differential expression
stats = data.frame()
for(s in sex){
  id = which(datMeta$Sex == s)
  meta = datMeta[id,]
  exp = datExp[,id]
  mod = model.matrix( ~Diagnosis + Age ,data = meta)
  
  fit <- lmFit(exp, mod)
  efit <- eBayes(fit, trend = T, robust = T)
  meta <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
  meta <- meta[,c(1,4,5)]
  colnames(meta) = c("logFC", "p.value", "fdr")
  sig = subset(meta, p.value<0.05) 
  sig$sex = s
  stats = rbind(stats,sig)
  
}

write.table(stats,"./results/tables/differential/CSF1_proteomics_sex_stats.txt",
            col.names = T,row.names = T)
}