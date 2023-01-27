#st2_networkParam.R
#---Here we compute the dendrogram and network parameters for the modules
rm(list=ls())

library(WGCNA); library(biomaRt);library(ggpubr); 
library(ComplexHeatmap); library(pheatmap);library(RColorBrewer)
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/bulk_analysis/replication_msbb/")
load("./codes/network/parameters/Combined_for_network.Rdata")

enableWGCNAThreads()
allowWGCNAThreads()
condition =TRUE

#load TOM comptued by rWGCNA
load("./codes/network/parameters/network_signed161_exprSet-block.1.RData")

geneTree = hclust(1-as.dist(TOM), method="average")

# Iterate WGCNA parameters for robustness -- this takes a while
colors = vector(mode="list")
labels = vector(mode="list")
if (condition){
    for (minModSize in c(20,50,100)) {
      for (dthresh in c(0.05,0.15, 0.25)) {
        for(ds in c(3:4)) { 
          print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
          
          tree = cutreeHybrid(
            dendro = geneTree,
            minClusterSize = minModSize,
            pamStage = FALSE, 
            cutHeight = 0.999,
            deepSplit = ds,
            distM = as.matrix(1-as.dist(TOM)))
          merged = mergeCloseModules(
            exprData = t(datExp), 
            colors = tree$labels, 
            cutHeight=dthresh)
          colors = cbind(colors, labels2colors(merged$colors))
          labels = c(labels, 
                     paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
        }
      }
    }
  }


save(file="./codes/network/parameters/WGCNA_diffParams.Rdata", geneTree, colors, labels)

for (i in 1:length(labels)){
  l = length(unique(colors[,i]))
  col = sort(as.character(unique(colors[,i])))
  lab = sort(labels2colors(0:(l-1)))
  if(all(col == lab)){
    print(c(i,l))
  }
  
}

if(condition){
  pdf("./results/figures/network/parameters/WGCNA_diffParams.pdf", width=18, height=14)
  
  par(mar= c(2,5,2,2))
  plotDendroAndColors(geneTree,
                      colors,guideHang = 0.2,
                      addGuide=T,
                      dendroLabels=F)
  
  plotDendroAndColors(geneTree, 
                      colors,guideHang = 0.2,
                      groupLabels = labels,
                      addGuide= TRUE,
                      dendroLabels=FALSE,
                      main="Dendrogram",
                      cex.colorLabels=0.5)
  dev.off()
}

## choosing parameters 
colors2 = cbind(colors[,16], colors)
labels = c("Final Modules", labels)

pdf("./results/figures/network/parameters/WGCNA_finalModules.pdf", width=10, height=14)
plotDendroAndColors(geneTree,
                    colors2,guideHang = 0.1,
                    groupLabels = labels,
                    addGuide=T,
                    setLayout = T,autoColorHeight = T,
                    dendroLabels=F,
                    cex.colorLabels=0.5)
dev.off()

# Finalized Parameters DS=4,MMS=100,DCOR=0.1,PAM=FALSE
wgcna_parameters = list(powers =  12)
wgcna_parameters$minModSize = 100
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 20000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 4  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = FALSE

if(condition){
  tree = cutreeHybrid(
    dendro = geneTree,
    minClusterSize = wgcna_parameters$minModSize,
    pamStage = wgcna_parameters$pamStage,
    cutHeight = 0.999, 
    deepSplit = wgcna_parameters$ds,
    distM = as.matrix(1-as.dist(TOM)))
  
  merged = mergeCloseModules(exprData= t(datExp),
                             colors = tree$labels,
                             cutHeight = wgcna_parameters$minHeight)
}

## change labels manually based on their order number
lab = sort(unique(merged$colors))
orderedLab = 0:(length(lab)-1)

for(i in 1:length(lab)){
  if(lab[i]!=orderedLab[i]){
    merged$colors[merged$colors==lab[i]] = orderedLab[i]
  }  
}

# collect colors 
colors = labels2colors(merged$colors)

# rename modules
mods = paste0("M",merged$colors)
n = length(unique(merged$colors))
names(mods) = colors
# obtain number of genes per module
colTable = as.data.frame(table(mods))

write.table(table(colors), 
            file = "./results/tables/network/ModuleColorsNumber.csv", row.names = F, col.names = T )

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
       x = "module")+
  theme_few(base_size = 10)+
  theme(legend.position = "none",
        text = element_text(size = 8, color = "black"),
        axis.text.x.bottom  =element_text(angle = 45, hjust = 1))
mod_p

ggsave("./results/figures/network/parameters/ModuleSize.pdf", 
       width = 5, height = 3,
       mod_p)

pdf(file = "./results/figures/network/parameters/finalModule.pdf", 
    height = 3, width = 8)

plotDendroAndColors(geneTree,
                    colors,
                    rowTextAlignment = "center", 
                    groupLabels = "Module",
                    main = "",
                    cex.colorLabels = 0.5,
                    addGuide = T,
                    dendroLabels = F)

dev.off()

MEs = moduleEigengenes(expr = t(datExp),
                       merged$colors, 
                       softPower = wgcna_parameters$powers,
                       scale = F)
kMEtable = signedKME(t(datExp),
                     MEs$eigengenes)
tableSup = data.frame(kMEtable)

######check correlation with metadata
MEs0 = as.data.frame(MEs$eigengenes)

# MEs0 = orderMEs(MEs0)
MEs0 = MEs0[,colnames(MEs0) != "ME0"]
colnames(MEs0) = gsub("ME","M",colnames(MEs0))

# merge datMeta and modules eigengene
ME.datMeta = cbind(datMeta, MEs0)

# a df of correlations
moduleCors = as.data.frame(matrix(NA, 
                                  nrow = ncol(datMeta), 
                                  ncol = ncol(MEs0)))
colnames(moduleCors) = colnames(MEs0);
rownames(moduleCors) = colnames(datMeta)

# exlcude sample and subjects from correlation 
moduleCors = moduleCors[c("Dx","Age","Sex","Apo1","Apo2","NP.1","CDR","PMI"),]

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

saveRDS("./codes/network/mod_cov_cor.rds",object = moduleCors)
saveRDS("./codes/network/mod_cov_P.rds",object = modulesPval)

#color
library(circlize)
h_col = colorRamp2(c(0,max(moduleCors)),c("white","blue"))
annot_mod = data.frame(module = colnames(moduleCors))
rownames(annot_mod)= colnames(moduleCors)
mod_col = unique(names(mods[mods != "M0"]))
names(mod_col) = unique(mods[mods !="M0"])
annot_col = list(module = mod_col)

## heatmap options
ht_opt( legend_border = "black",
        heatmap_border = TRUE,
        annotation_border = TRUE)

pdf("./results/figures/network/Modules_covariates_heatmap.pdf",
    width = 14, height = 7)

ComplexHeatmap::Heatmap(signif(moduleCors,2),
                        col = h_col,
                        top_annotation = HeatmapAnnotation(df = annot_mod,col = annot_col,show_legend = F),
                        name = "Correlation",
                        cluster_columns = T,
                        column_names_side = "top",
                        show_column_dend = F,
                        row_names_gp = gpar(fontsize = 8),
                        width = unit(30,"cm"),
                        height = unit(8,"cm"),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(modulesPval[i, j] < 0.05)
                            grid.text(signif(modulesPval[i, j],1), x, y,gp = gpar(fontsize = 7))
                        })
dev.off()


#----Annotating gene IDs
getinfo <- c( "ensembl_gene_id",
              "external_gene_name",
              "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
datProbes <- getBM(
  attributes = getinfo,
  filters = c("ensembl_gene_id"),
  values = rownames(tableSup),
  mart = mart)

tableSup = cbind(datProbes, 
                 data.frame(Module.Color=colors,
                            Module.name = paste0("M",merged$colors)),
                 tableSup)

write.csv(file="./results/tables/network/Modules_kMEtable.csv", tableSup)
save(file="./codes/network/parameters/finalNetwork.RData",
     datExp,
     datMeta,
     datProbes,
     geneTree,
     wgcna_parameters,
     colors,
     cols,mods,
     MEs,
     kMEtable)
