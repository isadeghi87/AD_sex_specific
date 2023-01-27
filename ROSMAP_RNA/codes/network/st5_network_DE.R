#st4_network_module_DE.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(WGCNA,ggplot2,reshape2,dplyr,magrittr,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
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
  
  ## filter significant ones
  id = unique(df$module[df$anova<0.05])
  df = df[df$module %in% id,]
  
  ## filter those correlated with cognition
  ## load correlation with covariates 
  load("./data/network/mod_cov_cor.Rdata")
  cordat = melt(moduleCors)
  colnames(cordat) = c("var","mod","cor")
  
  p = melt(modulesPval)
  cordat$p = p$value
  cordat$fdr = p.adjust(cordat$p,"fdr")
  
  cordat = subset(cordat,abs(cor)>0.2 & fdr <0.05)
  id = unique(as.character(cordat$mod))
  
  sexdat = subset(df,module %in% id )
  
  ## add sex-specific modules
  sexdiff = readRDS("./results/tables/network/sex_diff_modules.rds")
  sexdat$group = sexdiff$group[match(sexdat$module,sexdiff$mod)]
    
  sexdat$symbol = ""
  sexdat$symbol[sexdat$p<0.05] = "*"
  
  for(s in sexdat$group){
    dat = sexdat[sexdat$group==s,]
    
    col = dat$color
    names(col) = dat$module
    dat$module = factor(dat$module,levels = unique(dat$module))
    dat$diagnosis = factor(dat$diagnosis,levels = c("AsymAD","AD"))
    dat$symbol = ""
    dat$symbol[dat$fdr<0.05] = "*"
    
    
    ##  add anova value to the plot
    txt = dat[dat$anova<0.05 & dat$diagnosis =="AD",]
    txt$anova = signif(txt$anova,2)
    
    # plot
    p = ggplot(dat, aes(x = diagnosis,
                         y = beta,
                         group = diagnosis,
                         fill = module,
                         label = symbol)) + 
      facet_grid(sex~module,scales = "free") +
      geom_bar(stat = "identity",
               aes(fill = module),
               width = 0.5,show.legend = T,
               position = position_dodge(),
               color = "black") +
      geom_errorbar( aes(ymin = (beta - se), ymax = (beta + se)),
                     position = position_dodge(width = 0.5),
                     # size = 0.2,
                     width = 0.2) +
      geom_text(color = "red",
                size = 4,
                aes(y = beta + sign(beta) * se + sign(beta) * .01),
                position = position_dodge(width = 0.5)) +
      geom_text(data = txt,aes(label = paste0("P = ",anova),
                               y = 0.07, x= 1.5),size=4)+
      scale_fill_manual(values = col) +
      labs(y = "beta", x = "", 
           title = paste(s,"modules")) +
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
           filename = paste("./results/figures/network/sex_",s,"_cell_mod_beta.pdf",sep = ""),
           width = ifelse(s=="specific",6,14),
           height = ifelse(s=="specific",5,7))
    
  }
  
  saveRDS(sexdat,"./results/tables/network/sex_mod_DE.rds")
  
}
