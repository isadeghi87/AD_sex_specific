### tSNE.R
if(T){
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")

## libraries
pacman::p_load(ggplot2,Rtsne,ggthemes,ggpubr,limma,edgeR,
               ggalt,ggfortify,factoextra)

### load data
load("./codes/covariates/normalized_data.Rdata")

### PCA and MDS
set.seed(123)
sex = unique(datmeta$Sex)
plots = list()

for(s in sex){
  idx = which(datmeta$Sex == s)
  meta = datmeta[idx,]
  exp = datExp[,idx]
  
  # pc = prcomp(t(exp))
  # p = factoextra::fviz_pca_ind(pc,col.ind = meta$Dx,
  #                              # label="none"
  #                              )+
  #                            scale_color_brewer(palette="Set1")+ 
  #                            labs(title = s)+
  #                            theme_bw()+
  #                            theme(panel.grid = element_blank())
  # 
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
  theme_bw()+
  theme(panel.grid = element_blank())
p
plots[[s]] = p 
rm(p)
}

gg = ggarrange(plotlist = plots,common.legend = T,legend = "right")

ggsave("./results/figures/tsne/tsne.pdf",width = 6,height = 3,plot =gg)
}
