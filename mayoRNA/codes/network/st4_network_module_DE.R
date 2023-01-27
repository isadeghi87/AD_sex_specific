## network_module_DE.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  library(WGCNA);library(ggplot2); 
  library(reshape); library(nlme);library(dplyr)
  library(magrittr);library(ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/ROSMAP_RNA/")
  sex = c("female","male")
  datMeta$Diagnosis = gsub("AsymAD","Asym",datMeta$Diagnosis)
  datMeta$Diagnosis = factor(datMeta$Diagnosis,levels = c("CTL","Asym","AD"))
  
  load("./data/network/parameters/finalNetwork.RData")
  ## preserved modules
  idx = c("M1","M19")
  res = data.frame()
  
    
  
    eigmat = MEs$eigengenes
    colnames(eigmat) = gsub("ME","M",colnames(eigmat))
    n = ncol(eigmat)
    df = data.frame()
    
    for (s in sex){
      idx = which(datMeta$Sex==s)
      meta = datMeta[idx,]
      eig = eigmat[idx,]
      
      for (m in colnames(eigmat)) {
        me = as.numeric(eig[[m]])
        
        mod = lm(me ~ Diagnosis, data = meta)
        anov = anova(mod)["Diagnosis","Pr(>F)"]
        mod = as.data.frame(summary(mod)$coefficients)
        mod = mod[-1,-3]
        mod = cbind(mod,anov)
        colnames(mod) = c("beta","se","p","anova")
        mod$diagnosis = c("AsymAD","AD")
        mod$module = m
        mod$sex = s
        df = rbind(df,mod)
        
      }
    }
    
    df$fdr = p.adjust(df$p, "fdr")
    df$color = names(mods)[match(df$m,mods)]
    rownames(df) = NULL
    
    col = df$color
    names(col) = df$module  
    
    ## preserved module
    idx = "M19"
    dat = df[df$module==idx,]
    ##  add anova value to the plot
    txt = dat[dat$anova<0.05 & dat$diagnosis =="AD",]
    txt$anova = signif(txt$anova,2)
    
    # plot
    p = ggplot(dat, aes(x = diagnosis,
                        y = beta,
                        group = diagnosis,
                        fill = module)) + 
      facet_grid(sex~module,scales = "free") +
      geom_bar(stat = "identity",
               aes(fill = module),
               width = 0.3,show.legend = T,
               position = position_dodge(),
               color = "black") +
      geom_errorbar( aes(ymin = (beta - se), ymax = (beta + se)),
                     position = position_dodge(width = 0.5),
                     # size = 0.2,
                     width = 0.2) +
      geom_text(data = txt,aes(label = paste0("P = ",anova),
                               y = 0.07, x= 1.5),size=4)+
      scale_fill_manual(values = col) +
      labs(y = "beta", x = "", 
           title = "") +
      scale_x_discrete()+
      theme_bw(base_size = 14)+
      theme(legend.position = "none",
            strip.text = element_text(size = 14,face = "bold",color = "white"),
            axis.title = element_text(size = 12,angle = 90),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.background.y = element_rect(fill = "darkblue"),
            strip.background.x = element_rect(fill = "darkgreen"))
    
  p

  ggsave(plot = p,
         filename = "./results/figures/network/rosmap_rna_de.pdf",
         width = 4,
         height = 4)
}
