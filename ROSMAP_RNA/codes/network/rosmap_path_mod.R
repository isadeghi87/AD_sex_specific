## correlation of modules eigenvalue with ROSMAP pathological covariates
if(T){
rm(list=ls()); options(stringsAsFactors = F)
  
pacman::p_load(WGCNA,ggpubr,ComplexHeatmap,RColorBrewer,ggthemes)
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/ROSMAP_RNA/")
sex = c("female","male")

## load preserved modules
load("./data/network/parameters/finalNetwork.RData")
idx = "M19"
datMeta$cogn_global_lv = as.numeric(datMeta$cogn_global_lv)
datMeta$gpath = as.numeric(datMeta$gpath)

## filter ROSMAP
eig = MEs$eigengenes
colnames(eig) = gsub("E","",colnames(eig))
eig = eig[,idx]
all(rownames(eig)==rownames(datMeta))

## merge the data
merged = cbind(datMeta,eig)

res = data.frame(mod="NA",s="NA",var="NA",rho="NA",p="NA")
cov = c("gpath","cogn_global_lv")

for(mod in idx){
  for(s in sex){
    dat = merged[merged$Sex==s,]
    for(var in cov){
      x = dat$eig
      y = dat[,var,drop=T]
      corres = cor.test(x,y)
      rho = unname(corres$estimate)
      p = unname(corres$p.value)
      res = rbind(res,c(mod,s,var,rho,p))
    }
  }
}

res= res[-1,]
colnames(res) = c("mod","Sex","cov","rho","p")
res$p = as.numeric(res$p)
res$p = round(res$p,2)
res$rho = as.numeric(res$rho)
res$rho = round(res$rho,2)

fin.res = merge(merged,res,by = c("mod","Sex"))
fin.res$p[fin.res$Diagnosis != "AD"]= NA
fin.res$rho[fin.res$Diagnosis != "AD"]= NA

## add text
txt = fin.res[fin.res$Diagnosis =="AD",]
cols = c("gpath","cogn_global_lv")
fin.res = pivot_longer(merged,cols =cols ,
                       names_to = "vars",values_to = "cov_val")

## plot for global pathology and cognition
plots = list()
for(i in cols){
  dat = fin.res[fin.res$vars==i,]
  tx = res[res$cov==i,]
  pp = ggplot(dat, aes(x = cov_val, y = eig))+
    geom_point(aes(fill = Sex),shape = 21,color= "black")+
    geom_smooth(method = "lm",se=T,aes(color = Sex))+
    geom_text(data = tx,aes(label = paste0("P = ",p,"\n",
                                           "rho= ",rho),
                            color = Sex,group=Sex,
                            y = 0.15, x= 1.5),size=3,
              position = position_dodge(width = 5))+
    labs(x = i)+
    facet_grid(scales = "free")+
    theme_bw(base_size = 12)+
    theme(strip.text = element_text(size = 14,face = "bold",color = "white"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          # strip.background.y = element_rect(fill = "darkblue"),
          strip.background.x = element_rect(fill = "darkBlue"))
  print(pp)
  plots[[i]] = pp
  
}

gg = ggarrange(plotlist = plots,ncol = 1,
          labels = letters[1:length(plots)],
          common.legend = T,legend = "right")

ggsave("./results/figures/network/rosmap_rna_scatter_cor.pdf",
       plot =gg,width = 5,
       height = 7)

