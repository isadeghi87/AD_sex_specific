## overlap of genes between proteomics and transcriptome
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(gprofiler2,ggplot2,ggpubr,dplyr,ggthemes)
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  
  # Loading data ####
  load("./data/network/parameters/finalNetwork.RData")
  protein = tableSup
  
  ## load transcriptome network
  load("../ROSMAP_RNA/data/network/parameters/finalNetwork.RData")
  rna = tableSup
  
  ## overlap of M7 in protein and M19 in RNA
  protein = protein$gene[protein$mods=="M7"]
  rna = rna$gene[rna$mods=="M19"]
  
  ##  overlap
  ol.gene = intersect(protein,rna)
  write.table(ol.gene,"./results/tables/network/protein_RNA_intersect.txt")
  
  ## intersect with gene card AD gene list
  gc = read.csv("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/GeneCard_AD_genes.csv")
  idx = intersect(ol.gene,gc$Gene.Symbol)
  
  res = gc[gc$Gene.Symbol %in% idx,]
  res = res[,c(1:2)]
  rownames(res) = NULL
  write.table(res,"./results/tables/network/geneCard_AD_overlap.txt")
  
  sex = c("female","male")
  dat = datExp[ol.gene,]
  col = c("CTL"= "red","AsymAD"="blue","AD"="green")
  h.plots = list()
  
  for(s in sex){
    idx = which(datMeta$Sex ==s)
    meta = datMeta[idx,]
    exp = dat[,idx]
    
    #annotation of modules
    annot_mod = data.frame(Diagnosis = meta$Diagnosis)
    rownames(annot_mod) = colnames(exp)
    annot_col = list(Diagnosis = col)
    top_ann = HeatmapAnnotation(df = annot_mod,col = annot_col,
                                show_annotation_name = F)
    
    h = Heatmap(exp,
                top_annotation = top_ann,
                show_row_names = F,
                column_title = s,
                name = "expression",
                show_column_names = F)
    h.plots[[s]] = h
  }
  
  h.plots  

  
  }