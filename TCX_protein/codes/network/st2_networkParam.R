#st2_networkParam.R
#---Here we compute the dendrogram and network parameters for the modules
if(T){
rm(list=ls())

library(WGCNA); library(biomaRt);library(ggpubr); 
library(ComplexHeatmap); library(pheatmap);library(RColorBrewer)
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_temoral/")
load("./codes/network/parameters/Combined_for_network.Rdata")

# enableWGCNAThreads()
# allowWGCNAThreads()
load("./codes/tcx_proteomics_normalized.Rdata")
sex = unique(datMeta$Sex)

for(s in sex){
  # subsetting data for sex
  id = which(datMeta$Sex==s)
  exp = datExp[,id]
  meta = datMeta[id,]
  
  #load TOM comptued by rWGCNA
  network = list.files("./codes/network/",pattern = paste0("network_signed_",s,"_"))
  load(paste0("./codes/network/",network))
  
geneTree = hclust(1-as.dist(TOM), method="average")

# Iterate WGCNA parameters for robustness -- this takes a while
colors = vector(mode="list")
labels = vector(mode="list")

    for (minModSize in c(10,15,20)) {
      for (dthresh in c(0.01, 0.1,0.3)) {
        for(ds in c(3:4)) { 
         pam = FALSE
          print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
          
          tree = cutreeHybrid(
            dendro = geneTree,
            minClusterSize = minModSize,
            pamStage = pam, 
            minSplitHeight = dthresh,
            cutHeight = 1,
            deepSplit = ds,
            distM = as.matrix(1-as.dist(TOM)))
          # merged = mergeCloseModules(
          #   exprData = t(voomExp),
          #   colors = tree$labels,
          #   cutHeight=dthresh)
          colors = cbind(colors, labels2colors(tree$labels))
          labels = c(labels, 
                     paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
        }
      }
  }
  pdf(paste("./results/figures/network/WGCNA_diffParams_",s,".pdf",sep=""), width=18, height=14)
  
  par(mar= c(5,5,5,5))
  plotDendroAndColors(geneTree, 
                      colors,
                      guideHang = 0.2,
                      groupLabels = labels,
                      addGuide= TRUE,
                      dendroLabels=FALSE,
                      main="Dendrogram",
                      cex.colorLabels=0.5)
  dev.off()

## choosing parameters 
colors2 = cbind(colors[,2], colors)
labels = c("Final Modules", labels)

pdf(paste("./results/figures/network/WGCNA_finalModules_",s,".pdf",sep=""), width=10, height=14)
plotDendroAndColors(geneTree,
                    colors2,guideHang = 0.1,
                    groupLabels = labels,
                    addGuide=T,
                    setLayout = T,autoColorHeight = T,
                    dendroLabels=F,
                    cex.colorLabels=0.5)
dev.off()

# Finalized Parameters
power = stringr::str_extract(pattern = "[0-9]++",network)
power = as.numeric(power)
wgcna_parameters = list(powers =  power)
wgcna_parameters$minModSize = 10
wgcna_parameters$minHeight = 0.05
wgcna_parameters$bsize = 4000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 4  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = FALSE

  tree = cutreeHybrid(
    dendro = geneTree,
    minClusterSize = wgcna_parameters$minModSize,
    pamStage = wgcna_parameters$pamStage,
    minSplitHeight = 0.05,
    cutHeight = 1, 
    deepSplit = wgcna_parameters$ds,
    distM = as.matrix(1-as.dist(TOM)))
  
  # merged = mergeCloseModules(exprData= t(datExp),
                             # colors = tree$labels,
                             # cutHeight = wgcna_parameters$minHeight)


# collect colors 
colors = labels2colors(tree$labels)

# rename modules
mods = paste0("M",tree$labels)
n = length(unique(tree$colors))
names(mods) = colors

# obtain number of genes per module
colTable = as.data.frame(table(mods))

# colors
cols = table(names(mods)) %>% sort(decreasing = T) %>% names
colTable = colTable[order(colTable$Freq,decreasing = T),]
colTable$color = names(mods)[match(colTable$mods,mods)]
colTable$mods = factor(colTable$mods,levels = colTable$mods)
colTable$color = factor(colTable$color,levels = colTable$color)
library(ggthemes)

mod_p = ggplot(colTable, aes( x = mods,
                              y = Freq,
                              fill = mods))+
  geom_point(
    colour= "black",
    stroke=1,
    size = 4,shape = 21)+
  scale_fill_manual(values = levels(colTable$color))+
  scale_x_discrete()+
  scale_y_log10()+
  labs(y = "size",
       x = "module",
       title = paste(s,"modules"))+
  theme_few(base_size = 14)+
  theme(legend.position = "none",
        text = element_text(size = 8, color = "black"),
        axis.text.x.bottom  =element_text(angle = 45, hjust = 1))
mod_p

ggsave(paste("./results/figures/network/ModuleSize_",s,".pdf",sep=""), 
       width = 5, height = 3,
       mod_p)

## plot final module
pdf(file = paste("./results/figures/network/finalModule_",s,".pdf",sep=""), 
    height = 2, width = 7)

plotDendroAndColors(geneTree,
                    colors,
                    rowTextAlignment = "center", 
                    groupLabels = "Module",
                    main = "",
                    cex.colorLabels = 0.5,
                    addGuide = T,
                    dendroLabels = F)

dev.off()

MEs = moduleEigengenes(expr = t(exp),
                       tree$labels, 
                       softPower = wgcna_parameters$powers,
                       scale = F)
kMEtable = signedKME(t(exp),
                     MEs$eigengenes)
tableSup = data.frame(kMEtable)

######check correlation with metadata
MEs0 = as.data.frame(MEs$eigengenes)

# MEs0 = orderMEs(MEs0)
MEs0 = MEs0[,colnames(MEs0) != "ME0"]
colnames(MEs0) = gsub("ME","M",colnames(MEs0))

# merge datMeta and modules eigengene
ME.datMeta = cbind(meta, MEs0)

# a df of correlations
moduleCors = as.data.frame(matrix(NA, 
                                  nrow = ncol(meta), 
                                  ncol = ncol(MEs0)))
colnames(moduleCors) = colnames(MEs0);
rownames(moduleCors) = colnames(meta)


# exlcude sample and subjects from correlation 
moduleCors = moduleCors[c("Diagnosis","AgeAtDeath","ApoE.Risk",
                          "ApoE.Genotype","Braak"),]

#p values
modulesPval = moduleCors

#calculate linear regression between modules and covariates
for (mod in colnames(moduleCors)){
  for (var in rownames(moduleCors)){
    mcor =lm(ME.datMeta[,mod] ~ME.datMeta[,var], data = ME.datMeta )
    moduleCors[var,mod] = summary(mcor)$adj.r.squared
    modulesPval[var, mod] = anova(mcor)$`Pr(>F)`[1]
  }
}

###heatmap 
moduleCors = as.matrix(moduleCors)
modulesPval = as.matrix(modulesPval)
pval = signif(-log10(modulesPval),1)
symbol = moduleCors
symbol[modulesPval>0.05] =""
symbol[modulesPval<0.05] ="*"
symbol[modulesPval<0.01] ="**"
symbol[modulesPval<0.001] ="***"

saveRDS(paste("./codes/network/mod_cov_cor_",s,".rds",sep=""),object = moduleCors)
saveRDS(paste("./codes/network/mod_cov_P_",s,".rds",sep=""),object = modulesPval)


tableSup = cbind(data.frame(mods),gene,kMEtable)

write.csv(file="./results/tables/network/Modules_kMEtable.csv", tableSup)
save(file = paste("./codes/network/finalNetwork_",s,".RData",sep=""),
     exp,
     meta,
     gene,
     geneTree,
     wgcna_parameters,
     colors,
     cols,
     mods,
     MEs,
     kMEtable)
}
}
