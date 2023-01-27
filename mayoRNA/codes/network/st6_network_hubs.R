#st6_network_hubs.R
if(T){
  rm(list=ls())
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme);library(dplyr)
library(magrittr);library(ggthemes)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/")
sex = c("female","male")

# cell-type modules: 
res = data.frame()

for(s in sex){
#load final network
load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))

#### Make module eigengene-MDS plot ####
  kME = signedKME(t(exp), MEs$eigengenes)
  colnames(kME) = gsub("kME", "M", colnames(kME))
  id = which(colnames(kME)!="M0")
  kME = kME[,id]


maxsize = 10  #plot top 5 hub genes for each module

for(i in colnames(kME)){
  genes_idx = order(kME[which(mods == i),i], decreasing = T)[1:maxsize]
  hubs = rownames(exp[genes_idx,])
  symb = attr$gene_name[match(hubs,attr$gene_id)]
  res = rbind(res, cbind(s,i,symb))
}

}
write.csv(x = res,
          file = "./results/tables/network/sex_modul_hubGenes.csv")
##  hub genes for cell type mods #### 
# cell = cellmod$module[cellmod$sex==s]
# cons_kme = kME[,cell]
}
