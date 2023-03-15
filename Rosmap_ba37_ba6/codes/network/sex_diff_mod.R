## sex_diff_mod.R
#st4_network_module_DE.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(WGCNA,ggplot2,reshape2,dplyr,magrittr,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  sex = c("female","male")
  
  ## load preserved modules
  load("./data/network/parameters/finalNetwork.RData")
  
  eigmat = MEs$eigengenes
  colnames(eigmat) = gsub("ME","M",colnames(eigmat))
  datMeta$Sex = factor(datMeta$Sex,levels = c("female","male"))
  
  df = data.frame(m="NA",beta="NA",p="NA")
  for(m in unique(mods)){
    x = eigmat[,m]
    lm.test = summary(lm(x~Diagnosis*Sex,data = datMeta))
    p = lm.test$coefficients[2,4]
    beta =lm.test$coefficients[2,1] 
    df = rbind(df,c(m,beta,p))
  }
  
  df = df[-1,]
  df$p =as.numeric(df$p)
  df$beta =as.numeric(df$beta)
  colnames(df)[1] ="mod" 
  df$group = ifelse(df$p<0.05,"specific","shared")
  
  saveRDS(df,file = "./results/tables/network/sex_diff_modules.rds")
  
  ## colors for modules
  colname = unique(names(mods))
  names(colname)= unique(mods)
  
  ## plot
  pp = ggplot(df, aes(y = -log10(p),
                       x = reorder(mod,log10(p)),
                       fill = mod))+
    geom_bar(stat = "identity",colour = "black",
             width = 0.5 )+
    labs(x = "",
         title = "Modules with significant difference  between sexes")+
    scale_fill_manual(values =colname )+
    geom_hline(yintercept = -log10(0.05),lty=2)+
    theme_few(base_size = 12)+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,face = 'bold'),
          axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
          text = element_text(size = 10, color = "black"),
          strip.background = element_rect(fill = "black"),
          strip.text = element_text(color = "white",face = "bold",size=14))
  
  
  pp
    
  ggsave(plot = pp,
           filename = "./results/figures/network/sex_diff_modules.pdf",
           width = 7,
           height = 4)
    
  }
  
