### tSNE using top SMR genes
if(T){
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/")
  
  ## libraries
  library(pacman);
  p_load(ggplot2,factoextra,ggthemes,Rtsne)
  
  ### load data
  load("./codes/covariates/normalized_data.Rdata")
  
  ## read smr results
  smr = read_delim("../condition_overlap/disease_specific/AD_PA_CTL/data/gwas/eqtl_results/Brain_cis_eql_smr_results.txt") 
  smr$probeID = gsub("\\.[0-9]+","",smr$probeID)
  id = smr$probeID
  
  ### tSNE
  set.seed(123)
  sex = unique(datmeta$Sex)
  plots = list()
  
  for(s in sex){
    idx = which(datmeta$Sex == s)
    meta = datmeta[idx,]
    exp = datExp[,idx]
    index = intersect(id,rownames(exp))
    exp = exp[index,]
    
    ## PCA
    set.seed(123)
    tsne = Rtsne(X = t(exp),
                 # pca = T,
                 # pca_center  = T,
                 max_iter = 5000,perplexity = 20,theta = 0)
    res = as.data.frame(tsne$Y)
    rownames(res) = colnames(exp)
    res$condition = meta$Dx
    
    ### plot
    p = ggplot(data = res, aes(x  = V1, y = V2,
                               fill = condition))+
      geom_point(size = 3,shape = 21, color = "black")+
      labs(x = "tSNE1",y = "tSNE2",title = s)+
      scale_fill_manual(values = c("grey","blue","red"))+
      theme_bw()+
      theme(panel.grid = element_blank())
    p
    plots[[s]] = p 
    rm(p)
    
  }
  
  gg = ggarrange(plotlist = plots,common.legend = T,legend = "right")
  
  ggsave("./results/figures/tsne/tsne_top_smr.pdf",width = 6,height = 3,plot =gg )
}
