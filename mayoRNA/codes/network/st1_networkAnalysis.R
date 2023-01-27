#st1_networkAnalysis.R
## step1 WGCNA on mega analysis

if(T){
  rm(list=ls())
  library(WGCNA);library(limma);library(edgeR)
  library(stringr)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/")
  
  # load data
  load("./codes/covariates/normalized_data.Rdata")
  
  sex = unique(datmeta$Sex)
  
  for(s in sex){
    id = which(datmeta$Sex == s)
    exp = datExp[,id]
    meta = datmeta[id,]
    ## NETWORK ANALYSIS
    ## ----------------
    multiExpr = vector(mode="list", length=1)
    multiExpr[[1]] = list(data = as.data.frame(t(exp)))
    bsize = 20000
    nSets = 1
    powers = c(seq(1,9,by=1),seq(10,30,by=2))
    enableWGCNAThreads()
    allowWGCNAThreads()
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
    
    pdf(paste("./results/figures/network/WGCNA.softThresh",s,".pdf",sep = ""))
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
    
    wgcna_parameters = list(powers =  14)
    wgcna_parameters$minModSize = 30  
    wgcna_parameters$minHeight = 0.25
    wgcna_parameters$bsize = 20000  ##block size needs to be larger than dim(datExpr)[1]
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
                       power = wgcna_parameters$powers[n], 
                       mergeCutHeight= wgcna_parameters$minHeight, 
                       nThreads=4, 
                       saveTOMFileBase=paste("./codes/network/parameters/network_signed_",s, 
                                             wgcna_parameters, sep=""), 
                       saveTOMs=TRUE, 
                       minModuleSize= wgcna_parameters$minModSize,
                       pamStage=wgcna_parameters$pamStage, 
                       reassignThreshold=1e-6, 
                       verbose = 3, 
                       deepSplit=wgcna_parameters$ds)
  }
}

save(file = "./codes/network/parameters/Combined_for_network.Rdata",datmeta,datExp,attr)
