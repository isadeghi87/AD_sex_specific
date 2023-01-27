
## st8_AD_genes_enrichment.R
## enrichment of GWAS using MAGMA 
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(gplots,WGCNA,ggpubr,dplyr,magrittr,
                 GeneOverlap,tidyverse,ComplexHeatmap)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  # source("../condition_overlap/disease_specific/AD_PA_CTL/bulk_analysis/codes/source/Fisher_overlap.R")
  
  ## read magma results
  magma = read.table("C:/Users/crgcomu/Desktop/Iman/Brain_meta/tools/magma/results/snp.magma_output.genes.out",header = T) 
  magma = subset(magma,P<0.05) ## keep significant
  attr = read_rds("C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/gencode28_attributes.rds")
  magma$symbol = attr$gene_name[match(magma$GENE,attr$gene_id)]
  id= unique((magma$symbol))
  
  
  load("./data/network/parameters/finalNetwork.RData")
  modules = unique(mods)
  
  res = data.frame()
  #Calculate spearman's correlation between gene module membership and GWAS gene significance
  
    # for(m in modules) {
    #   col = paste("kME", m, sep="")
    #   col = gsub("kMEM","kME",col)
    #   id = which(tableSup$mods==m)
    #   set = intersect(magma$symbol,gene[id])
    #   
    #   x  = -log10(magma$P[match(set, magma$symbol)])
    #   y = kMEtable[match(set,tableSup$gene), col]
    #   cor = cor.test(x,y,method="spearman")
    #   r = as.numeric(cor$estimate)
    #   p = cor$p.value
    #   res = rbind(res,cbind(m,r,p))
    # }
    # res$p = as.numeric(res$p)
    # res$r = as.numeric(res$r)
    # res$fdr = p.adjust(res$p,"fdr")
    # res$fdr = as.numeric(res$fdr)
    # res$color = names(mods)[match(res$m,mods)]
    # res$m = factor(res$m,levels = res$m)
    # res$symbol = ""
    # res$symbol[res$fdr<0.05]="*"
    # res$symbol[res$fdr<0.01]="**"
    # res$symbol[res$fdr<0.001]="***"
  
  
  ## use over-representation
  for(m in modules) {
    id1 = gene[mods==m]
    id2 = id # genes from magma
    n = length(union(id1,id2))
    o.obj = GeneOverlap::newGeneOverlap(id1,id2)
    tst = testGeneOverlap(o.obj)
    # odd ratio
    or = as.numeric(signif(getOddsRatio(tst),2))
    # get P value
    p= as.numeric(signif(getPval(tst),2))
    res = rbind(res,cbind(m,or,p))
  }
  
  res$p = as.numeric(res$p)
  res$or = as.numeric(res$or)
  res$fdr = p.adjust(res$p,"fdr")
  res$color = names(mods)[match(res$m,mods)]
  # res$m = factor(res$m,levels = res$m)
  res$symbol = ""
  res$symbol[res$fdr<0.05]="*"
  res$symbol[res$fdr<0.01]="**"
  res$symbol[res$fdr<0.001]="***"
  
  res$log.fdr = -log10(res$fdr)
  
  ## colors
  colManuel = names(mods)
  names(colManuel) = mods
  
  ##plot
  p   = ggplot(res, aes(x = reorder(m,-log10(p)),
                        y = log.fdr,
                        fill = m,
                        label = symbol))+
    geom_bar(stat = "identity",
             colour = "black",width = 0.5 )+
    # facet_wrap(~s,ncol=1,scales = "free")+
    geom_abline(intercept = -log10(0.05),slope=0,lty=2)+
    # geom_text(aes(label = symbol, y = -log10(fdr) +0.1),
    #           color = "red",size = 5, show.legend = F)+
    scale_fill_manual(values = colManuel)+
    labs(x = "",
         title= "AD GWAS enrichment",
         y = "-log10(fdr)")+
    theme_bw(base_size = 10)+
    theme(axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
          strip.background = element_rect(fill = "darkblue"),
          plot.title = element_text(hjust = 0.5,face = "bold"),
          strip.text = element_text(size = 8,face = "bold",color = "white"))+
    guides(fill = "none", size = "none")
  p
  ggsave("./results/figures/network/AD_gwas_magma_enrichment.pdf",
         width = 7,height = 4,plot = p)
  
}  
