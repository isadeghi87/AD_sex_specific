#st1_networkAnalysis.R
## step1 WGCNA on mega analysis for TMT consensus

if(T){
  rm(list=ls())
pacman::p_load(WGCNA,limma,edgeR,biomaRt,stringr)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/ROSMAP_RNA/")

# load data
load("./data/normalized.Rdata")
gene = row.names(datExp)
sex = unique(datMeta$Sex)

## NETWORK ANALYSIS
## ----------------
multiExpr = vector(mode="list", length=1)
multiExpr[[1]] = list(data = as.data.frame(t(datExp)))
bsize = 16000
nSets = 1
powers = c(seq(1,9,by=1),seq(10,30,by=2))
# enableWGCNAThreads()
# allowWGCNAThreads()
n=1

## Compute soft-threshold for mega-analysis
multiExpr[[n]]$softThresh = pickSoftThreshold(data= multiExpr[[n]]$data, 
                                              networkType = "signed", 
                                              corFnc = "bicor",
                                              verbose = 5,
                                              powerVector = powers,
                                              blockSize = bsize)

sft = multiExpr[[n]]$softThresh
sft$powerEstimate
power = sft$powerEstimate
max.sft = max(sft$fitIndices$SFT.R.sq)
# power = sft$fitIndices$Power[sft$fitIndices$SFT.R.sq==max.sft]

pdf("./results/figures/network/WGCNA.softThresh.pdf")
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

## DS=4,MMS=15,DCOR=0.05
wgcna_parameters = list(powers =  power)
wgcna_parameters$minModSize = 15  
wgcna_parameters$minHeight = 0.2
wgcna_parameters$bsize = 16000  ##block size needs to be larger than dim(datExpr)[1]
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
                   minKMEtoStay = 0.3,
                   power = wgcna_parameters$powers[n], 
                   mergeCutHeight= 0.99,
                   minSplitHeight = 0.05,
                   nThreads=4, 
                   saveTOMFileBase=paste("./data/network/parameters/network_signed", 
                                         power, sep="_"), 
                   saveTOMs=TRUE, 
                   minModuleSize= wgcna_parameters$minModSize,
                   pamStage=wgcna_parameters$pamStage, 
                   reassignThreshold=1e-6, 
                   verbose = 3, 
                   deepSplit=wgcna_parameters$ds)
save(file = "./data/network/parameters/Combined_for_network.Rdata",datMeta,datExp,gene)
}

