
## st8_AD_genes_enrichment.R
## enrichment of GWAS using MAGMA 
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(gplots,WGCNA,ggpubr,dplyr,magrittr,
                 GeneOverlap,tidyverse,ComplexHeatmap)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific//")
  source("../condition_overlap/disease_specific/AD_PA_CTL/bulk_analysis/codes/source/Fisher_overlap.R")
  
  ## read magma results
  magma = read.table("C:/Users/crgcomu/Desktop/Iman/Brain_meta/tools/magma/results/snp.magma_output.genes.out",header = T) 
  magma = subset(magma,P<0.05) ## keep significants
  id= unique((magma$GENE))
  
  # cell-type modules: 
  cellmod = readRDS("./codes/network/cell_modules.rds")
  sex = c("female","male")
  enrich = data.frame()
  
  for(s in sex){
    #load final network
    load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
    modules = unique(mods[mods!="M0"])
    
    res = data.frame()
    #Calculate spearman's correlation between gene module membership and GWAS gene significance
    
    #   for(m in modules) {
    #     col = paste("kME", m, sep="")
    #     col = gsub("kMEM","kME",col)
    #     genes = intersect(magma$GENE, rownames(exp))
    #     x  = -log10(magma$P[match(genes, magma$GENE)])
    #     y = kMEtable[match(genes,rownames(kMEtable)), col]
    #     cor = cor.test(x,y,method="spearman")
    #     r = as.numeric(cor$estimate)
    #     p = cor$p.value
    #     res = rbind(res,cbind(s,m,r,p))
    #   }
    #   res$p = as.numeric(res$p)
    #   res$r = as.numeric(res$r)
    #   res$fdr = p.adjust(res$p,"fdr")
    #   res$fdr = as.numeric(res$fdr)
    #   res$color = names(mods)[match(res$m,mods)]
    #   res$m = factor(res$m,levels = res$m)
    #   res$symbol = ""
    #   res$symbol[res$fdr<0.05]="*"
    #   res$symbol[res$fdr<0.01]="**"
    #   res$symbol[res$fdr<0.001]="***"
    #   enrich = rbind(enrich,res)
    # }
    
    
    ## use over-representation
    for(m in modules) {
      id1 = rownames(exp)[mods==m]
      id2 = id # genes from magma
      n = length(union(id1,id2))
      o.obj = GeneOverlap::newGeneOverlap(id1,id2)
      tst = testGeneOverlap(o.obj)
      # odd ratio
      or = as.numeric(signif(getOddsRatio(tst),2))
      # get P value
      p= as.numeric(signif(getPval(tst),2))
      res = rbind(res,cbind(s,m,or,p))
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
    enrich = rbind(enrich,res)
  }
  
  enrich$log.fdr = -log10(enrich$fdr)
  
  ## colors
  colManuel = names(mods)
  names(colManuel) = mods
  
  ##plot
  p   = ggplot(enrich, aes(x = reorder(m,-log.fdr),
                           y = log.fdr,
                           fill = m,
                           label = symbol))+
    geom_bar(stat = "identity",
             colour = "black",width = 0.5 )+
    facet_wrap(~s,ncol=1,scales = "free")+
    geom_abline(intercept = -log10(0.05),slope=0,lty=2)+
    # geom_text(aes(label = symbol, y = -log10(fdr) +0.1),
    #           color = "red",size = 5, show.legend = F)+
    scale_fill_manual(values = colManuel)+
    labs(x = "",
         y = "-log10(fdr)")+
    theme_bw(base_size = 10)+
    theme(axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
          strip.background = element_rect(fill = "darkblue"),
          strip.text = element_text(size = 8,face = "bold",color = "white"))+
    guides(fill = "none", size = "none")
  p
  ggsave("./results/figures/network/AD_gwas_magma_enrichment.pdf",
         width = 6,height = 4,plot = p)
  
}  
