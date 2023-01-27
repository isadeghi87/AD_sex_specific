## st7_proteomics_enrichment.R
rm(list=ls()); options(stringsAsFactors = F)
library(gplots); library(biomaRt); library(WGCNA)
library(ggpubr); library(pheatmap); library(dplyr); library(magrittr)
library(GeneOverlap); library(factoextra);library(FactoMineR)
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

## annotate ensemble ID for datexp
getinfo <- c("ensembl_gene_id",
             "external_gene_name")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
datexp.attr <- getBM(
  attributes = getinfo,
  filters = "ensembl_gene_id",
  values = rownames(datExp),
  mart = mart)
datexp.attr = datexp.attr[match(rownames(datExp),datexp.attr$ensembl_gene_id),]

### proteomics modules from Johnson et al , Nat Med.2020
protein_dat = readxl::read_xlsx("./data/proteomics_Erik_et_al/supplementary.xlsx",
                        sheet = 15, skip =2) 
protein_dat$UniqueID = str_split(string = protein_dat$UniqueID, pattern = "\\|",simplify = T)[,1]

## get ensemble id for proteins
getinfo <- c("ensembl_gene_id",
             "external_gene_name")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl") ## Gencode v28
bm <- getBM(
  attributes = getinfo,
  filters = "external_gene_name",
  values = protein_dat$UniqueID,
  mart = mart)

bm = bm[match(protein_dat$UniqueID,bm$external_gene_name),]
protein_dat$ensID = bm$ensembl_gene_id[match(protein_dat$UniqueID,bm$external_gene_name)]

# remove grey module
protein_dat = protein_dat[protein_dat$`net$colors`!="grey",]

# rename protein modules
modN = 1:13
labelcol = labels2colors(1:13)
names(labelcol)= paste0("pM",modN)
protein_dat$module = names(labelcol)[match(protein_dat$`net$colors`,labelcol)]

## check overlap for proteomics modules ####
protein_overlap = matrix(NA, nrow = length(unique(protein_dat$module)),
                      ncol = ncol(kME))
rownames(protein_overlap) = unique(protein_dat$module)
colnames(protein_overlap) = colnames(kME)
protein_overlap.or = protein_overlap


for (m in colnames(protein_overlap)) {
  for(p in rownames(protein_overlap)) {
    f = ORA(rownames(datExp)[mods==m],
           protein_dat$ensID[protein_dat$module==p],
            rownames(datExp),bm$ensembl_gene_id )
    protein_overlap.or[p,m] = as.numeric(f[[1]])
    protein_overlap[p,m] = as.numeric(f[[2]])
  }
}

protein_overlap.fdr = p.adjust(protein_overlap,"fdr")
dim(protein_overlap.fdr) = dim(protein_overlap); 
dimnames(protein_overlap.fdr) = dimnames(protein_overlap)
protein_overlap.fdr[protein_overlap.or<1] = 1

if(condition){
## protein annotation
protein_annot = data.frame(proteomics = rownames(protein_overlap))
rownames(protein_annot) = rownames(protein_overlap)
protein_annot$proteomics= str_replace_all(protein_annot$proteomics,
                                          c("pM1" = "Neuron",
                                            "pM4" = "Astro/Microglia",
                                            "pM2" ="Oligo",
                                            "pM5"= "Microglia/Endo",
                                            "pM7"= "Endo"))
p_cellMod= c("pM1","pM2","pM4","pM5","pM7")
protein_annot$proteomics[!rownames(protein_annot) %in% p_cellMod ] = ""

## annotation for our modules
mod_ann  = data.frame(transcriptomics = colnames(protein_overlap))
rownames(mod_ann) = colnames(protein_overlap)
p_cellMod= c("M1","M4","M6","M9","M11","M13","M15","M5","M16")
mod_ann$transcriptomics[!rownames(mod_ann) %in% p_cellMod ] = ""
mod_ann$transcriptomics = str_replace_all(mod_ann$transcriptomics,
                                c(
                                  "M4" = "Microglia",
                                  "M5" = "Neuron",
                                  "M6"= "Endo",
                                  "M9" = "Neuron",
                                  "M11" = "Neuron",
                                  "M13" = "Astro",
                                  "M15" = "Astro",
                                  "M16" = "Neuron"))
mod_ann$transcriptomics[mod_ann$transcriptomics=="M1"] = "Oligo"




## annotation colors
ann_col = list(transcriptomics = c( Neuron = "magenta",
                                    Astro = "midnightblue",
                                    Oligo = "turquoise",
                                    Microglia = "yellow",
                                    Endo = "red"),
               proteomics = c( Neuron = "magenta",
                               `Astro/Microglia` = "midnightblue",
                               Oligo = "turquoise",
                               `Microglia/Endo` = "yellow",
                               Endo = "red"))
## heatmap
pdf("./results/wgcna/proteomics_module_overlap.pdf",width = 10,height = 4)
ComplexHeatmap::pheatmap(protein_overlap.or,
                         color = redblue(1000)[500:1000],
                         annotation_colors  = ann_col,
                         treeheight_col = 4,treeheight_row = 4,
                         annotation_row = protein_annot,
                         annotation_names_row = T,angle_col = '45',
                         annotation_names_col = T,
                         cellwidth = 25,cellheight = 15,
                         annotation_col = mod_ann,
                         name = "overlap ratio",
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           if(protein_overlap.fdr[i, j] < 0.05)
                             grid.text(signif(protein_overlap.fdr[i, j],1), x, y,gp = gpar(fontsize = 7))
                         })
dev.off()
}


# Enrichment for inflammatory genes for astro and microgia ####
astroMarkers = as.data.frame(read_xlsx("./data/proteomics_Erik_et_al/supplementary.xlsx",
                                       sheet = 22, skip = 2))
astroP = matrix(NA, nrow = ncol(astroMarkers),
                                  ncol = ncol(kME))
rownames(astroP) = colnames(astroMarkers)
colnames(astroP) = colnames(kME)
astro.or = astroP

for (c in colnames(astroP)) {
  for(r in rownames(astroP)) {
    id1 = datexp.attr$external_gene_name[mods ==c]
    idx = !is.na(unique(astroMarkers[,r]))
    id2 = unique(astroMarkers[,r])[idx]
    o.obj = GeneOverlap::newGeneOverlap(id1,id2)
    tst = testGeneOverlap(o.obj)
    # odd ratio
    astro.or[r,c] = signif(getOddsRatio(tst),2)
    # get P value
    astroP[r,c] = signif(getPval(tst),2)
  }
}

astro.fdr = p.adjust(astroP,"fdr")
dim(astro.fdr) = dim(astroP)
colnames(astro.fdr) = colnames(astroP)
rownames(astro.fdr) = rownames(astroP)
astro.fdr[astro.or<1] = 1

astroLabels = c("A1\ninflammatory","A2\ntissue repair","A1/A2\nshared")

## heatmap

## check microglia markers ####
micMarkers = as.data.frame(read_xlsx("./data/proteomics_Erik_et_al/supplementary.xlsx",
                                       sheet = 25,
                                     skip = 2))
colnames(micMarkers) = c("Anti-inflammatory","Pro-inflammatory","Homeostatic")
markers = union(micMarkers$`Anti-inflammatory`, micMarkers$`Pro-inflammatory`)
markers = union(markers,micMarkers$Homeostatic)

## annotate 
getinfo <- c("ensembl_gene_id",
             "external_gene_name","hgnc_symbol")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl") ## Gencode v28
## df for results
micP = matrix(NA, nrow = ncol(micMarkers),
                ncol = ncol(kME))
rownames(micP) = colnames(micMarkers)
colnames(micP) = colnames(kME)
mic.or = micP

for (c in colnames(micP)) {
  for(r in rownames(micP)) {
    id1 = rownames(datExp)[mods ==c]
    idx = !is.na(unique(micMarkers[,r]))
    idx = unique(micMarkers[,r])[idx]
    
    micAttr <- getBM(
      attributes = getinfo,
      filters = "external_gene_name",
      values = idx,
      mart = mart)
    id2 = micAttr$ensembl_gene_id
    o.obj = GeneOverlap::newGeneOverlap(id1,id2,
                                        spec = "mm9.gene")
    tst = testGeneOverlap(o.obj)
    # odd ratio
    mic.or[r,c] = signif(getOddsRatio(tst),2)
    # get P value
    micP[r,c] = signif(getPval(tst),2)
  }
}

mic.fdr = p.adjust(micP,"fdr")
dim(mic.fdr) = dim(micP)
colnames(mic.fdr) = colnames(micP)
rownames(mic.fdr) = rownames(micP)
mic.fdr[mic.or<1] = 1

micLabels = c("Anti-inflammatory","Pro-inflammatory","Homeostatic")


### draw heatmaps ####
colModule = cols[-1]
names(colModule) = mods[match(cols[-1],names(mods))]

annot_col = data.frame(module = colnames(astroP))
rownames(annot_col)= colnames(astroP)
annot_color = list(module = colModule)
top_ann = HeatmapAnnotation(df=annot_col,
                            col =annot_color,show_legend = F,show_annotation_name = F)
if(condition){
## heatmaps
h1 = Heatmap(-log10(astro.fdr),
             col = rev(redblue(1000))[500:900],
             row_names_side = "left",row_labels = astroLabels,
             top_annotation = top_ann,
             show_column_names = T,column_names_side = "top",column_names_rot = 45,
             column_names_gp = gpar(fontsize = 8),
             column_dend_side = "top",
             row_names_gp = gpar(fontsize = 8),
             width = unit(18,"cm"),
             height = unit(4,"cm"),
             # row_title = "Astrocyte",row_title_side = "right",
             border = "black",
             cluster_columns = T,
             show_row_dend = F,
             name = "astrocyte\n-log10(FDR)",
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(astro.fdr[i, j] < 0.05)
                 grid.text(signif(astro.fdr[i, j],1), x, y,gp = gpar(fontsize = 8))
             })

## heatmap of correlation
h2 = Heatmap(-log10(mic.fdr),
             col = redblue(1000)[500:900],
             row_names_side = "left",row_labels = micLabels,
             show_column_names = F,
             row_names_gp = gpar(fontsize = 8),
             width = unit(18,"cm"),
             height = unit(3,"cm"),
             # row_title = 4Microglia",row_title_side = "right",
             border = "black",
             show_column_dend = F,show_row_dend = F,
             name = "microglia\n-log10(FDR)",
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(mic.fdr[i, j] < 0.05)
                 grid.text(signif(mic.fdr[i, j],1), x, y,gp = gpar(fontsize = 8))
             })




## heamap annotation of modules

pdf("./results/wgcna/astro_micro_enrichment.pdf", width = 9,height = 4)
hlist =  h1 %v% h2 
draw(hlist)
dev.off()
}
