#st6_network_hubs.R
if(T){
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme);library(dplyr)
library(magrittr);library(ggthemes)
condition = TRUE

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
sex = c("female","male")
# cell-type modules: 
cellmod = readRDS("./codes/network/cell_modules.rds")
res = data.frame()

for(s in sex){
#load final network
load(paste0("./codes/network/finalNetwork_",s,".RData",sep=""))

#### Make module eigengene-MDS plot ####
  eigmat = MEs$eigengenes
  colnames(eigmat) = gsub("ME","M",colnames(eigmat))
  eigmat = eigmat[,-1]
  kME = signedKME(t(exp), MEs$eigengenes)
  colnames(kME) = gsub("kME", "M", colnames(kME))
  kME = kME[,-1]

##  hub genes for cell type mods #### 
cell = cellmod$module[cellmod$sex==s]
cons_kme = kME[,cell]

maxsize = 5  #plot top 5 hub genes for each module
hub.list = list()

for(i in colnames(cons_kme)){
  genes_idx = order(cons_kme[which(mods == i),i], decreasing = T)[1:maxsize]
  hubs = rownames(exp[genes_idx,])
  symb = datProbes$gene_name[match(hubs,datProbes$gene_id)]
  hub.list[[i]] = symb
}

hubGenes = do.call("cbind",hub.list)
hubGenes = melt(hubGenes)
hubGenes$sex = s
hubGenes = hubGenes[,2:4]
colnames(hubGenes)[1:2]= c("module","gene")
res = rbind(res,hubGenes)

}
}
write.csv(x = res,
          file = "./results/tables/network/sex_modul_hubGenes.csv")

