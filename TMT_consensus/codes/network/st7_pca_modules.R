
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
sex = c("female","male")

## load preserved modules
load("./data/network/parameters/finalNetwork.RData")

sex = ("female","male")
set = c("M45","M53")

  for(s in sex){
    idx = which(datMeta$Sex == s)
    exp = datExp[,idx]
    metadata = datMeta[idx,]
  
  for (m in set){
    res = exp[id,]
    pc = prcomp(t(dat),scale. = T,center = F)
    tsn = Rtsne::Rtsne(pc$x,pca = F,set.seet = 123)
    tsnRes = as.data.frame(tsn$Y)
    rownames(tsnRes) = colnames(dat)
    tsnRes$Sex = s
    tsnRes$module = m
    
    tsnRes$diagnosis = metadata$Diagnosis
    res = rbind(res,tsnRes)
    
  }
  }


p = ggplot(res,aes(x = V1,y = V2))+
  geom_point(aes(color = diagnosis),size = 2)+
  facet_grid(Sex~module)+
  labs(x = "tSNE1",y = "tSNE2",title = "")+
  scale_color_brewer(palette="Set1")+ 
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold",size = 14))

ggsave(filename = "./results/figures/network/pca_mods.pdf",plot = p,
       width = 8,height = 5)  
