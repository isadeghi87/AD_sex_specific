###--MSBB.R
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/bulk_analysis/replication_msbb/")

#libraries
library(reshape); library(stringr); library(dplyr);library(edgeR);
library(readr);  library(tidyverse);library(stringr);library(sva);
library(WGCNA);library(limma);library(nlme)

#options
options(stringsAsFactors = F)

# Load covariates meta Data
datMeta = read.csv("./metadata/MSBB_Covariates_analyzed.csv")
rownames(datMeta) = datMeta$Sample
datMeta$Diagnosis = gsub("Normal","CTL",datMeta$Diagnosis)
datMeta$Dx = factor(datMeta$Diagnosis,levels = c("CTL","AD"))
datMeta$Sex = as.factor(datMeta$Sex)
datMeta$Race = as.factor(datMeta$Race)
datMeta$batch = as.factor(datMeta$batch)
datMeta$region = as.factor(datMeta$Brain_Region)
datMeta$lobe = gsub("Prefrontal","Frontal",datMeta$Brain_Lobe)
datMeta$lobe = as.factor(datMeta$lobe)
datMeta$RIN[is.na(datMeta$RIN)] = mean(datMeta$RIN, na.rm=T)
datMeta$Study = "MSBB"
datMeta$Subject_ID = datMeta$individualIdentifier

##----Load Expression Data
datExpr <- read.csv("./expression_files/MSBB_ExprData.csv", sep = "\t",
                        row.names = 1)
rownames(datExpr) = gsub("\\.\\d{1,2}", "", rownames(datExpr))

idx = match(rownames(datMeta),colnames(datExpr))
idx = na.omit(idx)
datExpr = datExpr[,idx]
idx = match(colnames(datExpr),rownames(datMeta))
idx = na.omit(idx)
datMeta = datMeta[idx,]

# remove duplicates
id = which(colnames(datExpr) %in% c("hB_RNA_7765","hB_RNA_8215","hB_RNA_10742"))
datExpr = datExpr[,-id]
datMeta = datMeta[-id,]

## keep only Temporal lobe
id = which(datMeta$lobe == "Temporal")
datExpr = datExpr[,id]
datMeta = datMeta[id,]
datMeta$lobe = factor(datMeta$lobe, levels = unique(datMeta$lobe))
datMeta$region = factor(datMeta$region, levels = unique(datMeta$region))

## Remove outliers based on network connectivity z-scores
normadj <- (0.5+0.5*bicor(datExpr))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,col=c("blue" , "green"),pch=19, xlab = "",ylab= "Z-score")
legend("bottomleft",legend = levels(datMeta$Dx), 
       col = c("blue", "green"),pch=19,cex=.7)
abline(h=-2, lty=2)
outliers = (z.ku < -2)
table(outliers)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

all(colnames(datExpr)==rownames(datMeta))

#### filter lowly expressed genes ####
id = rowSums(cpm(datExpr)>0.5)> 0.3*ncol(datExpr)
datExpr = datExpr[id,]

# normalize counts and fitting
design = model.matrix(~Dx  , data = datMeta)
colnames(design) = gsub("Dx","", colnames(design))

#### DE ####
datExpr = datExpr+1
v = voom(datExpr, design = design)
datExpr =v$E 


## remove batch
datMeta$batch = factor(datMeta$batch,levels = unique(datMeta$batch))
datExpr = removeBatchEffect(datExpr,batch = datMeta$batch,design = design)

## tSNE
library(Rtsne);library(ggplot2);library(ggthemes)
set.seed(123)
tsn = Rtsne(t(datExpr),theta = 0,pca = T,verbose = 3)
dat = as.data.frame(tsn$Y)
rownames(dat) = colnames(datExpr)
dat$condition = datMeta$Dx
dat$region = datMeta$region
dat$batch = datMeta$batch
  
ggplot(dat,aes(x = V1,y = V2, color = condition))+
  geom_point()

fit = lmFit(v, design = design)
efit <- eBayes(fit, trend = T, robust = T)

#### Alzheimer ####
ad.stats <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
ad.stats <- ad.stats[,c(1,4,5)]
colnames(ad.stats) = c("logFC", "p.value", "fdr")

write.csv(x = ad.stats, file = "./results/tables/MSBB_AD_stats.csv")
save(file="./codes/MSBB_normalized.Rdata",datMeta,datExpr)
