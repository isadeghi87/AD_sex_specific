#st_network_hubs.R
if(T){
  rm(list=ls())
  pacman::p_load(WGCNA,ggplot2,reshape,dplyr,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  sex = c("female","male")
  
  # cell-type modules: 
  
    #load final network
    load("./data/network/parameters/finalNetwork.RData")
    
    #### Make module eigengene-MDS plot ####
    kME = signedKME(t(datExp), MEs$eigengenes)
    colnames(kME) = gsub("kME", "M", colnames(kME))
    
    n = 10  #plot top 5 hub genes for each module
  res = data.frame()
    for(i in colnames(kME)){
      genes_idx = order(kME[which(mods == i),i], decreasing = T)[1:n]
      hubs = rownames(datExp[genes_idx,])
      symb = gene[match(hubs,rownames(datExp))]
      res = rbind(res, cbind(i,symb))
    }
    
  write.csv(x = res,
            file = "./results/tables/network/TMT_hubGenes.csv")
}
