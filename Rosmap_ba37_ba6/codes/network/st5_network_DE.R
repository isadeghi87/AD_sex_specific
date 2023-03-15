#st4_network_module_DE.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(WGCNA,ggplot2,reshape2,dplyr,magrittr,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/Rosmap_ba37_ba6/")
  sex = c("female","male")
  
  ## load preserved modules
  load("./data/network/parameters/finalNetwork.RData")
  
  ## load network 
  eigmat = MEs$eigengenes
  colnames(eigmat) = gsub("ME","M",colnames(eigmat))
  
  #Step 1 - Calculate module-trait P values, logFC, and SEM
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
  
  ## preserved matched modules
  pres = read.csv("../TMT_consensus/results/tables/network/preservation/rosmap_BA37.csv")
  pres = pres[pres$z>1.96,]
  
  ## specific modules from TMT
  spec = readRDS("../TMT_consensus/results/tables/network/sex_mod_DE.rds")
  spec = spec[spec$group=="specific",]
  id = unique(spec$module)
  
  ## sex specific that are preserved 
  idx = intersect(pres$labels,id)
  matchmod = read.csv("../TMT_consensus/results/tables/network/preservation/rosmap_BA37_matched_modules.csv")
  matchmod = matchmod[matchmod$TMT==idx,]  
  
  ### M7 is M9 in tcx
  ## and M39 is M13 
  df = df[df$module %in% c("M9","M13"),]
  df$tmt_mod = matchmod$TMT[match(df$module,matchmod$Tcx)]
  
  ## add colors of TMT network
  load("../TMT_consensus/data/network/parameters/finalNetwork.RData")
  df$tmt_color = names(mods)[match(df$tmt_mod,mods)]
  
  col = df$tmt_color
  names(col) = df$tmt_mod
  
  # dat$tmt_mod = factor(dat$tmt_mod,levels = unique(dat$tmt_mod))
  df$diagnosis = factor(df$diagnosis,levels = c("AsymAD","AD"))
  df$symbol = ""
  df$symbol[df$fdr<0.05] = "*"
  
  ##  add anova value to the plot
  txt = df[df$anova<0.05 & df$diagnosis =="AD",]
  txt$anova = signif(txt$anova,2)
  
  # plot
  p = ggplot(df, aes(x = diagnosis,
                     y = beta,
                     group = diagnosis,
                     fill = tmt_mod,
                     label = symbol)) + 
    facet_grid(sex~tmt_mod,scales = "free") +
    geom_bar(stat = "identity",
             aes(fill = tmt_mod),
             width = 0.5,show.legend = T,
             position = position_dodge(),
             color = "black") +
    geom_errorbar( aes(ymin = (beta - se), ymax = (beta + se)),
                   position = position_dodge(width = 0.5),
                   # size = 0.2,
                   width = 0.2) +
    # geom_text(color = "red",
    #           size = 4,
    #           aes(y = beta + sign(beta) * se + sign(beta) * .01),
    #           position = position_dodge(width = 0.5)) +
    geom_text(data = txt,aes(label = paste0("P = ",anova),
                             y = 0.07, x= 1.5),size=4)+
    scale_fill_manual(values = col) +
    labs(y = "beta", x = "",
         title = "BA37 temporal cortex proteomics") +
    scale_x_discrete()+
    theme_bw(base_size = 14)+
    theme(legend.position = "none",
          strip.text = element_text(size = 14,face = "bold",color = "white"),
          axis.title = element_text(size = 12,angle = 90),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          plot.title = element_text(size = 12,hjust = 0.5, face = "bold"),
          strip.background.y = element_rect(fill = "darkblue"),
          strip.background.x = element_rect(fill = "darkgreen"))
  
  p
  ggsave(plot = p,
         filename = paste("./results/figures/network/sex_mod_DE.pdf",sep = ""),
         width = 3,height = 5)
  
}


