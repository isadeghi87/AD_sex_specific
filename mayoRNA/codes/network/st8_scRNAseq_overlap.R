## st8_scRNAseq_overlap.R
## data from scRNA-seq of Grubman et al Nature neuro
## https://doi.org/10.1038/s41593-019-0539-4

rm(list=ls()); options(stringsAsFactors = F)
library(gplots); library(biomaRt); library(WGCNA)
library(ggpubr); library(dplyr); library(magrittr)
library(GeneOverlap); library(ComplexHeatmap)
library(readxl);library(tidyverse)

condition = TRUE
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/")
source("./codes/Fisher_overlap.R")

#---load final network
load("./codes/wgcna/finalNetwork.RData")

eigmat = MEs$eigengenes
eigmat = eigmat[,-1]
colnames(eigmat) = gsub("ME","M",colnames(eigmat))
kME = signedKME(t(datExp), MEs$eigengenes[,-1])
colnames(kME) = gsub("kME", "M", colnames(kME))

# saveRDS("./tables/attributes.rds",object = datexp.attr)
datexp.attr = read_rds("./tables/attributes.rds")

### scRNAseq results markers
sc_dat = readxl::read_xlsx("../data/scRNAseq/scRNAseq_Grubman.xlsm",
                                sheet = 2, skip =5) 
sc_dat$group = gsub("mg","microglia",sc_dat$group)

## check overlap for  modules ####
resP = matrix(NA, nrow = length(unique(sc_dat$group)),
                         ncol = ncol(kME))
rownames(resP) = unique(sc_dat$group)
colnames(resP) = colnames(kME)
res.or = resP


for (c in colnames(resP)) {
  for(r in rownames(resP)) {
    id1 = datexp.attr$external_gene_name[mods ==c]
    idx = subset(sc_dat,group == r)
    id2 = unique(idx$geneID)
    o.obj = GeneOverlap::newGeneOverlap(id1,id2, spec = "hg19.gene")
    tst = testGeneOverlap(o.obj)
    # odd ratio
    res.or[r,c] = as.numeric(signif(getOddsRatio(tst),2))
    # get P value
    resP[r,c] = as.numeric(signif(getPval(tst),2))
  }
}

res.fdr = p.adjust(resP,"fdr")
dim(res.fdr) = dim(resP)
colnames(res.fdr) = colnames(resP)
rownames(res.fdr) = rownames(resP)
res.fdr[res.or<1] = 1


colModule = cols[-1]
names(colModule) = mods[match(cols[-1],names(mods))]

annot_col = data.frame(module = colnames(resP))
rownames(annot_col)= colnames(resP)
annot_color = list(module = colModule)
top_ann = HeatmapAnnotation(df=annot_col,
                            col =annot_color,show_legend = F,show_annotation_name = F)
# if(condition){
  ## heatmaps
  h1 = Heatmap(-log10(res.fdr),
               col = rev(redblue(1000))[500:900],
               row_names_side = "left",
               top_annotation = top_ann,
               show_column_names = T,
               column_names_side = "top",column_names_rot = 45,
               column_names_gp = gpar(fontsize = 8),
               column_dend_side = "top",
               row_names_gp = gpar(fontsize = 8),
               width = unit(18,"cm"),
               height = unit(4,"cm"),
               border = "black",
               cluster_columns = T,
               show_row_dend = F,
               name = "-log10(FDR)",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(resP[i, j] < 0.05)
                   grid.text(signif(res.fdr[i, j],1), x, y,gp = gpar(fontsize = 8))
               })
  

pdf("./results/wgcna/scRNA_new_marker_overlap.pdf",
    width = 8.5,height = 3)
draw(h1)
dev.off()


sc = c_dat = readxl::read_xlsx("./data/scRNAseq/scRNAseq_Grubman.xlsm",
                               sheet = 4, skip =5)
sc = subset(sc,cond =="AD")
