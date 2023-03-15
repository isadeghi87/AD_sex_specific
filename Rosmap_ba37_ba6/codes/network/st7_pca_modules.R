## pca_modules.R
# here we compute tSNE for sex-specific module
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
pacman::p_load(ggplot2,Rtsne,ggthemes,umap)

## load preserved modules
load("./data/network/parameters/finalNetwork.RData")

sex = ("female","male")
set = c("M7")
idx = which(mods=="M7")
dat = datExp[idx,]

if(T){
res = data.frame()
  
for(s in sex){
    idx = which(datMeta$Sex == s)
    exp = dat[,idx]
    metadata = datMeta[idx,]
  
    # tsn = Rtsne::Rtsne(t(exp),pca = T,
    #                    max_iteration = 2000,
    #                    initial_dims=20,
    #                    theta = 0,
    #                    perplexity = 20,
    #                    normalize = F,
    #                    pca_center = F,
    #                    set.seet = 123,)
    umdat = umap(t(exp),preserve.seed = T)
    umRes = as.data.frame(umdat$layout,)
    # rownames(tsnRes) = colnames(exp)
    umRes$Sex = s
    
    umRes$diagnosis = metadata$Diagnosis
    res = rbind(res,umRes)
    
  }


p = ggplot(res,aes(x = V1,y = V2))+
  geom_point(aes(color = diagnosis),size = 2)+
  facet_wrap(~Sex,scale = "free")+
  labs(x = "UMAP 1",y = "UMAP 2",
       title = "")+
  scale_color_brewer(palette="Set1")+ 
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold",size = 14))

p
}
ggsave(filename = "./results/figures/network/pca_mods.pdf",plot = p,
       width = 6,height = 3)  
