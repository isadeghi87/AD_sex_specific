## gene enrichment analysis for each module
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  library(pacman)
  p_load(gprofiler2,ggplot2,ggpubr,dplyr,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  attr = readRDS("C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/gencode28_attributes.rds")
  
    # Loading data ####
    load("./data/network/parameters/finalNetwork.RData")
    
    ## filter each module
  res = data.frame()
    for (m in unique(mods)) {
      id = which(mods == m)
      idx = gene[id]
      query = attr$gene_id[match(idx,attr$gene_name)]
      query = query[!is.na(query)]
      ## Calculate GO enrichment of Disease-associated modules
      go = gost(query,
                organism = "hsapiens",
                significant = FALSE,
                # user_threshold = 0.05,
                correction_method = "fdr",
                domain_scope = "annotated",
                sources = c("KEGG","GO:BP","GO:MF"))
      
      # go = gProfileR::gprofiler(query,significant = F,
      #                           correction_method = "fdr",
      #                           src_filter  = c("KEGG","GO:BP"))

      go = go$result
      go = cbind(m,go)
      ## order results by p value
      go = go[order(go$p_value),]
      go = go[1:5,]
      res = rbind(res, go)
      
  }
  
  mod_col = unique(names(mods))
  names(mod_col) = unique(mods)
  
  ## sex specific modules
  sexmod = readRDS("./results/tables/network/sex_mod_DE.rds")
  idx = unique(sexmod$module[sexmod$group=="specific"])
  
  sexSpec = res[res$m %in% idx, ]

 p = ggplot(sexSpec,aes(x = -log10(p_value),y = term_name,fill = m))+
    geom_bar(stat = "identity",show.legend = F)+
    facet_wrap(~m,ncol = 1,scales = "free_y")+
    labs(y="")+
    theme_few()
 
 ggsave("./results/figures/network/go_sex_specific.pdf",p,
          width = 6,height = 4)
  
}
