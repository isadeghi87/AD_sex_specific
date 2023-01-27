#st1_networkAnalysis.R
## step1 WGCNA on mega analysis

rm(list=ls())
library(WGCNA);library(limma);library(edgeR);library(biomaRt)
library(stringr)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/bulk_analysis/replication_proteomics_CSF1/")

# load data
load("./codes/CSF1_proteomics_normalized.Rdata")

## NETWORK ANALYSIS
## ----------------
multiExpr = vector(mode="list", length=1)
multiExpr[[1]] = list(data = t(datExp))
bsize = 1000
nSets = 1
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()
n=1

## check genes for missing values 
id = goodGenes(datExpr = t(datExp),minFraction = 0.5,)
datExp = datExp[id,]

## Compute soft-threshold for mega-analysis
multiExpr[[n]]$softThresh = pickSoftThreshold(data= multiExpr[[n]]$data, 
                                              networkType = "signed", 
                                              corFnc = "bicor",
                                              verbose = 5,
                                              powerVector = powers,
                                              blockSize = bsize)

sft = multiExpr[[n]]$softThresh
sft$powerEstimate

if(T){
  pdf("./results/figures/network/parameters/WGCNA.softThresh.pdf")
  par(mfrow=c(1,2))
  n = 1
  plot(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       xlab="Soft Threshold Power",
       ylab=expression("Scale free R "^2),   
       type="n")
  text(sft$fitIndices[,1],
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers,
       cex = 0.7,
       col="blue",  
       xlab="Soft Threshold Power",
       ylab="Scale free R^2")
  abline(h= 0.85, col="black")
  plot(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       xlab = "Soft threshold power",
       ylab = "Mean connectivity",
       type = "n")
  text(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       labels = powers,
       cex = 0.7, 
       col="black")
  dev.off()
}

## power= 18
## DS=4,MMS=10,DCOR=0.05
wgcna_parameters = list(powers =  18)
wgcna_parameters$minModSize = 10  
wgcna_parameters$minHeight = 0.05
wgcna_parameters$bsize = 1000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 4  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = TRUE

n = 1

## ROBUST WGCNA 
## --> we want to generate a "consensus TOM" based on resampling datExp so as to make 
## --> module definitions robust to outliers
multiExpr[[n]]$netData = 
  blockwiseModules(datExpr=multiExpr[[n]]$data, 
                   maxBlockSize=wgcna_parameters$bsize,
                   networkType=wgcna_parameters$networkType, 
                   corType = wgcna_parameters$corFnc,
                   # minKMEtoStay = 0.3,
                   power = wgcna_parameters$powers[n], 
                   mergeCutHeight= 1,
                   minSplitHeight = wgcna_parameters$minHeight,
                   nThreads=4, 
                   saveTOMFileBase=paste("./codes/network/parameters/network_signed", 
                                         wgcna_parameters,as.character(n), "_exprSet", sep=""), 
                   saveTOMs=TRUE, 
                   minModuleSize= wgcna_parameters$minModSize,
                   pamStage=wgcna_parameters$pamStage, 
                   reassignThreshold=1e-6, 
                   verbose = 5,
                   deepSplit=wgcna_parameters$ds)


save(file = "./codes/network/parameters/Combined_for_network.Rdata",datMeta,datExp,attr)
