## csf1_prot_covariates.R
## covariates of csf1  proteomics
if(T){
  rm(list = ls())
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF1/")
  
  ## libraries
  library(pacman);
  p_load(ggplot2,factoextra,ggthemes,Rtsne,ggpubr)
  
  ## load data ####
  load("./codes/CSF1_proteomics_normalized.Rdata")
  
  
  ## plot covariates ####
  cols = c("grey","green")
  sex = unique(datMeta$Sex)
  plots = list()
  
  for(s in sex){
    idx = which(datMeta$Sex == s)
    meta = datMeta[idx,]
    exp = datExp[,idx]
    
    ## PCA
    pc = prcomp(t(exp))
    
    p = factoextra::fviz_pca_ind(pc,
                                 label="none",
                                 col.ind  = meta$Diagnosis)+
      labs(title = s)+
      scale_color_brewer(palette="Set1")+ 
      theme_bw()+
      theme(panel.grid = element_blank())
    
    plots[[s]] = p 
    
  }
  
  gg = ggarrange(plotlist = plots,common.legend = T,legend = "right")
  
  ggsave("./results/figures/pca/csf1_sex_pca.pdf",width = 7,height = 3,plot =gg )
  
}
