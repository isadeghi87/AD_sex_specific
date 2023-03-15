
## enrichment of BrainMeta eqtl and mendelian randomization (smr)
## results

if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  library(gplots)
  library(dplyr);library(ggrepel)
  library(GeneOverlap);library(ggthemes)
  library(tidyverse)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  source("../bulk_analysis/codes/source/Fisher_overlap.R")
  
  smr = read_delim("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/gwas/eqtl_results/Brain_cis_eql_smr_results.txt") 
  smr$probeID = gsub("\\.[0-9]+","",smr$probeID)
  sex = c("female","male")
  enrich = topsnp = txt= data.frame()
  
    #load final network
    load("./codes/network/parameters/finalNetwork.RData")
    # modules = unique(mods[mods!="M0"])
    df = data.frame(Gene = gene,
                    module = mods,
                    color = names(mods))
    df = merge(df,smr,by = "Gene")
    
    colManuel = names(mods)
    names(colManuel) = mods
    colManuel = colManuel[colManuel %in% df$color]
    
    p = ggplot(df, aes(x = b_SMR,
                       y = -log10(p_SMR_multi)))+
      geom_point(aes(fill = module),shape = 21, 
                 colour = "black",size = 3)+
      # geom_hline(yintercept = -log10(0.05),lty=2)+
      geom_text_repel(aes(label = Gene),size = 2.5)+
      scale_fill_manual(values = colManuel)+
      labs(x = "SMR beta",
           title = "Top significant gene from mendelian randomization",
           # subtitle = paste(n,"significant SNPs in total"),
           y = "-log10(SMR fdr)")+
      theme_bw(base_size = 10)+
      theme(axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
            strip.background = element_rect(fill = "#004D40"),
            strip.text = element_text(size = 12,face = "bold",color = "white"))
    
    ggsave("./results/figures/network/smr_enrichment.pdf",p,
           width = 6,height = 4)
    
    # df = subset(df,module !="M0")
    topsnp = rbind(topsnp,df)
    ## data for labels
    df = df[order(df$p_SMR_multi,decreasing = F),]
    txt = rbind(txt,df[1:10,])
    
    res = data.frame()
    #Calculate enrichment between gene module membership and GWAS gene significance
    for(m in modules) {
      id1 = gene[mods==m]
      id2 = smr$Gene
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
    enrich = rbind(enrich,res)
  }
  
  enrich$log.fdr = -log10(enrich$fdr) * sign(enrich$or) 
  ## only positive enrichment
  enrich$log.fdr[enrich$log.fdr<0]=0
  
  ## colors
  n = length(smr$topSNP)
  
  ##plot
  p = ggplot(enrich, aes(x = reorder(m,-log.fdr),
                           y = log.fdr,
                           fill = m))+
    geom_bar(stat = "identity",
             colour = "black",width = 0.5 )+
    facet_wrap(~s,ncol=1, scales = "free_x")+
    geom_abline(intercept = -log10(0.05),slope=0,lty=2)+
    scale_fill_manual(values = colManuel)+
    labs(x = "",
         title = "Enrichment for mendelian randomization of brain eQTLs",
         subtitle = paste(n,"significant SNPs in total"),
         y = "-log10(fdr)")+
    theme_bw(base_size = 10)+
    theme(axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
          strip.background = element_rect(fill = "#004D40"),
          strip.text = element_text(size = 10,face = "bold",color = "white"))+
    guides(fill = "none", size = "none")
  p
  ggsave("./results/figures/network/AD_smr_enrichment.pdf",
         width = 6,height = 5,plot = p)
  
  
  ## plot for top SNPs
  topsnp$module = as.character(topsnp$module)
  col = topsnp$color;names(col)=topsnp$module
  
   pp = ggplot(topsnp, aes(x = b_SMR,
                           y = -log10(p_SMR_multi)))+
     geom_point(aes(fill = module),shape = 21, 
                colour = "black",size = 3)+
     facet_wrap(~sex, ncol=1, scales = "free")+
     geom_hline(yintercept = -log10(0.05),lty=2)+
     geom_text_repel(data =txt ,aes(label = Gene),size = 2.5)+
     scale_fill_manual(values = col)+
     labs(x = "SMR beta",
          title = "Top significant gene from mendelian randomization",
          # subtitle = paste(n,"significant SNPs in total"),
          y = "-log10(SMR fdr)")+
     theme_bw(base_size = 10)+
     theme(axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
           strip.background = element_rect(fill = "#004D40"),
           strip.text = element_text(size = 12,face = "bold",color = "white"))
   
   ggsave("./results/figures/network/AD_smr_topSNP.pdf",
          width = 7, height = 5,plot = pp)
   
   write.csv(topsnp,"./results/tables/network/topSNP_smr.csv")

