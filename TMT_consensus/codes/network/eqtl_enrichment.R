
## eqtl_enrichment.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(gplots,WGCNA,ggpubr,dplyr,magrittr,
                 GeneOverlap,tidyverse,ComplexHeatmap)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  # source("../condition_overlap/disease_specific/AD_PA_CTL/bulk_analysis/codes/source/Fisher_overlap.R")
  
  ## read ROsMAP eQTL
  eqtl = read_delim("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/proteomics/ROSMAP/raw_data/ROSMAP_eQTLs.txt") 
  eqtl = subset(gwas,pValue<5e-8)
  id= unique(eqtl$featureName)
  
  ## load network data
  load("./data/network/parameters/finalNetwork.RData")
  modules = unique(mods)
  
  ## use over-representation
  res = data.frame()
  for(m in modules) {
    id1 = gene[mods==m]
    id2 = id # genes from magma
    (int = intersect(id1,id2))
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
  p   = ggplot(res, aes(x = reorder(m,-log.fdr),
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
         title= "eQTL enrichment",
         y = "-log10(fdr)")+
    theme_bw(base_size = 10)+
    theme(axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
          strip.background = element_rect(fill = "darkblue"),
          plot.title = element_text(hjust = 0.5,face = "bold"),
          strip.text = element_text(size = 8,face = "bold",color = "white"))+
    guides(fill = "none", size = "none")
  p
  ggsave("./results/figures/network/eqtl_enrichment.pdf",
         width = 7,height = 4,plot = p)
  
}  
