#st7_geneEnrichment.R

rm(list=ls()); options(stringsAsFactors = F)
library(gprofiler2); library(ggplot2); library(biomaRt); library(WGCNA)
library(ggpubr); library(ComplexHeatmap); library(dplyr); library(magrittr)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
sex = c("female","male")

# cell-type modules: 
cellmod = list(female = c("M1","M15","M2","M17","M18","M37","M39","M22","M9","M26","M12"),
               male = c("M1","M13","M3","M6","M17","M24","M9","M2","M20"))

plot.list = list()
for(s in sex){
  #load final network
  load(paste0("./codes/network/finalNetwork_",s,".RData",sep=""))
  
  cell = cellmod[[s]]
  all_mod = unique(mods)
  
  ## Calculate GO enrichment of Disease-associated modules
  res = data.frame()
  plotData = data.frame()
  for (module in all_mod) {
    # GO Enrichment
    query = rownames(exp)[mods == module]
    go = gost(query,
               organism = "hsapiens",
               ordered_query = FALSE,
               significant = F,
               user_threshold = 0.05,
               correction_method = "fdr",
               domain_scope = "annotated",
               sources = c("KEGG","GO","REAC","HP"))
    
    go = as.data.frame(go$result)
    go = cbind(module,go)
    ## order results by p value
    go = go[order(go$p_value),]
    go = go[1:4,]
    res = rbind(res, go)
  }
  
  mod_col = unique(names(mods))
  names(mod_col) = unique(mods)
  
  ## plot go for all modules
  p = ggplot(res, aes(x = -log(p_value),
                      y = reorder(term_name,log10(p_value)),
                      fill = module))+
    geom_bar(stat = "identity",show.legend = F)+
    facet_wrap(~module,ncol = 4,scales = "free")+
    scale_fill_manual(values = mod_col)+
    labs(x = "-log10(FDR)",y = "")+
    theme_bw(base_size = 9)
  
  ggsave(filename = paste0("./results/figures/network/all_modules_go_",s,".pdf"),
                           plot = p,
         width = 15, height = 15)
  
  #  plot data for cell type modules
  plotData = res[res$module %in% cell, c("module","term_name","p_value")]
  modCol = names(mods)[match(cell,mods)]
  names(modCol) = cell
  
  ## bar plot 
  p = ggplot(plotData, aes(x = -log10(p_value), 
                           y = reorder(term_name,log10(p_value)), 
                           fill = module))+
    geom_bar(stat = "identity",width = 0.4, color = "black",show.legend = F)+
    facet_wrap(~module, ncol = 1, scales = "free")+
    scale_fill_manual(values = modCol)+
    scale_x_log10()+
    labs(x = "-log10(FDR)", y = "", title = s)+
    theme_few(base_size = 12)+
    theme(strip.background.x = element_rect(fill = "#363062"),
          strip.text.x = element_text(color = "white",face = "bold",size = 12))+
    guides(y = guide_axis(position = "right"))
  
  plot.list[[s]] = p
}

pp = ggarrange(plotlist = plot.list,ncol = 2)

ggsave(pp,filename = "./results/figures/network/cell_mod_go.pdf",
         width = 15,height = 12)
  