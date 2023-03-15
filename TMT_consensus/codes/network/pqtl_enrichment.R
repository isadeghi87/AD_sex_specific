
##pqtl_enrichment.R
## enrichment for ROSMAP pQTL 
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(gplots,WGCNA,ggpubr,dplyr,magrittr,
                 GeneOverlap,tidyverse,ComplexHeatmap,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  
  ## read magma results
  pqtl = read.table("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/proteomics/ROSMAP/raw_data/ROSMAP.for_smr.pQTLv2.txt",header = T) 
  pqtl = subset(pqtl,P < 5e-8) ## keep significants
  pqtl$symbol = str_split(pqtl$Protein_GeneSymbol,pattern = "\\.",n=2,simplify = T)[,2]
  id= unique(pqtl$symbol)
  
  load("./data/network/parameters/finalNetwork.RData")
  modules = unique(mods)
  
  res = data.frame()
  
  ## use over-representation
  for(m in modules) {
    id1 = gene[mods==m]
    id2 = id 
    int = intersect(id1,id2)
    
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
  res$color = names(mods)[match(res$m,mods)]
  
  ## sex specific modules
  sexdif = read_rds("./results/tables/network/sex_mod_DE.rds")
  id = unique(sexdif$module)
  res = res[res$m %in% id, ]

  ## colors
  colManuel = names(mods)
  names(colManuel) = mods
  
  ##plot
  p   = ggplot(res, aes(x = reorder(m,-or),
                        y = or,
                        fill = m))+
    geom_bar(stat = "identity",
             colour = "black",width = 0.5 )+
    # facet_wrap(~s,ncol=1,scales = "free")+
    geom_abline(intercept = 1.96,slope=0,lty=2)+
    # geom_text(aes(label = symbol, y = -log10(fdr) +0.1),
    #           color = "red",size = 5, show.legend = F)+
    scale_fill_manual(values = colManuel)+
    labs(x = "",
         title= "AD pQTL enrichment",
         y = "Enrichment")+
    theme_bw(base_size = 14)+
    theme(axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
          strip.background = element_rect(fill = "darkblue"),
          plot.title = element_text(hjust = 0.5,face = "bold"),
          strip.text = element_text(size = 10,face = "bold",color = "white"))+
    guides(fill = "none", size = "none")
  p
  ggsave("./results/figures/network/pqtl_enrichment.pdf",
         width = 7,height = 4,plot = p)
  
  ## top SNPs for specific modules
  specMod = unique(sexdif$module[sexdif$group=="specific"])
  snp.res = manh.res = data.frame()
  
  for(m in specMod ){
    id = gene[mods==m]
    vars = pqtl[pqtl$symbol %in% id,]
    vars = vars[order(abs(vars$Beta),decreasing = F),]
    manh.res = rbind(manh.res,vars)
    # vars = vars[1,]
    # snp.res = rbind(snp.res,vars)
    
  }
  
  ## filter spec modules
  attr.mod = tableSup[tableSup$mods %in% specMod,]
  manh.res$module = attr.mod$mods[match(manh.res$symbol,attr.mod$gene)]
  
  ## manual colors
  colname = unique(names(mods[mods %in% specMod ]))
  names(colname) = unique(mods[mods %in% specMod ])
  
  # library(ggrepel)
  # pp = ggplot(manh.res,aes(x = Beta,y = -log10(P), fill = module))+
  #   geom_point(shape = 21,color = "black",size =3)+
  #   scale_fill_manual(values = colname)+
  #   geom_text_repel(data = snp.res,aes(label = SNP ))+
  #   theme_few()
  # 
  # ggsave("./results/figures/network/topSNPs_mod.pdf",plot = pp,
  #        width = 7,height = 4)
 write.table(manh.res,"./results/tables/network/topQTLs_mod.txt",sep = "\t")
  
  ## for M7 
  dat = manh.res[manh.res$module =="M7",]
  dat = dat[order(dat$Beta,decreasing = T),]
  idx = c("SNP","Chr","A1","A2","symbol","Beta","P","module")
  dat = dat[1:5,idx]
  rownames(dat) = NULL
  p = ggtexttable(dat,theme = ttheme(base_style = "mBlue"))
  ggsave("./results/figures/network/topSNPs_table.pdf",plot = p,
         width = 7,height = 4)
  
  ## genes that are overlap between protein and RNA
  ol = read.table("./results/tables/network/geneCard_AD_overlap.txt")
  ol.qtl = manh.res[manh.res$symbol %in% ol$Gene.Symbol,]
  write.table(ol.qtl)
  ## gwas
  # gwas = readxl::read_xlsx("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/gwas/AD_gwas_genes.xlsx",
  #                          sheet = "Supplementary Table 8a",skip = 5) 
  # gwas = gwas[gwas$`P-VALUE` < 5e-8,]## keep significants
  # colnames(gwas)[1] = "symbol"
  # ol.gwas= gwas[gwas$symbol %in% ol$Gene.Symbol, ]
 
  
}  
