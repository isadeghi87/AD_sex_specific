## mod_covariate_cor.R
## we compute correlation between covariates 
## and module eigenvalue

if(T) {
  rm(list=ls())
  pacman::p_load(WGCNA,ggpubr,ComplexHeatmap,RColorBrewer,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  
  load("./data/network/parameters/finalNetwork.RData")
  
  # datMeta$Braak =as.factor(datMeta$Braak)
  # datMeta$CERAD =as.factor(datMeta$CERAD)
  
  ######check correlation with metadata
  eigmat = as.data.frame(MEs$eigengenes)
  
  # eigmat = orderMEs(eigmat)
  colnames(eigmat) = gsub("ME","M",colnames(eigmat))
  
  #keep only ROSMAP
  id = which(!is.na(datMeta$study.ROSMAP))
  meta = datMeta[id,]
  eig = eigmat[id,]
  
  # merge datMeta and modules eigengene
  df = cbind(meta, eig)
  df$Sex = as.numeric(as.factor(df$Sex))
  
  # a df of correlations
  moduleCors = as.data.frame(matrix(NA, 
                                    nrow = ncol(datMeta), 
                                    ncol = ncol(eigmat)))
  colnames(moduleCors) = colnames(eigmat);
  rownames(moduleCors) = colnames(datMeta)
  
  # covariates
  moduleCors = moduleCors[c("Braak","CERAD","gpath","cogn_global_lv",
                            "amyloid","nft","tangles"),]
  
  #p values
  modulesPval = moduleCors
  
  #calculate linear regression between modules and covariates
  for (mod in colnames(moduleCors)){
    for (var in rownames(moduleCors)){
      mcor =cor.test(x = df[,mod],
                     y = df[,var], data = df )
      moduleCors[var,mod] = mcor$estimate
      modulesPval[var, mod] = mcor$p.value
    }
  }
  
  moduleCors = as.matrix(moduleCors)
  modulesPval = as.matrix(modulesPval)
  
  save(file = paste("./data/network/mod_cov_cor.Rdata",sep=""),moduleCors,modulesPval)
}