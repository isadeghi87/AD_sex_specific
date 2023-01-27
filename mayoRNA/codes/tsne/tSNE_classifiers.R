### tSNE_classifiers.R
if(T){
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  
  ## libraries
  library(pacman);
  p_load(ggplot2,factoextra,ggthemes,Rtsne)

    ### load data
  load("./codes/covariates/normalized_data.Rdata")
  ## read top important genes
  gene = read.csv("./results/tables/conditions_top_important_genes.csv")
  
  ### tSNE
  set.seed(123)
  sex = unique(datmeta$Sex)
  plots = list()
  

for(s in sex){
    idx = which(datmeta$Sex == s)
    meta = datmeta[idx,]
    exp = datExp[,idx]
    id = gene$var[gene$s==s]
    exp = exp[id,]
    
    ## PCA
    pc = prcomp(t(exp))
    p = factoextra::fviz_pca_ind(pc,label="none",col.ind  = meta$Dx)+
      labs(title = s)+
      scale_color_brewer(palette="Set1")+ 
      theme_bw()+
      theme(panel.grid = element_blank())
  
    plots[[s]] = p 
    
  }
  
  gg = ggarrange(plotlist = plots,common.legend = T,legend = "right")
  
  ggsave("./results/figures/tsne/tsne_classifiers.pdf",width = 7,height = 3,plot =gg )
}
