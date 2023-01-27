#st4_network_module_DE.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(WGCNA,ggplot2,reshape2,dplyr,magrittr,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/consensus/")
  sex = c("female","male")
  
  ## load preserved modules
  load("./data/network/parameters/finalNetwork.RData")
  res = data.frame()
  
  
  ## load network 
  eigmat = MEs$eigengenes
  colnames(eigmat) = gsub("ME","M",colnames(eigmat))
  id = which(colnames(eigmat)!="M0")
  eigmat = eigmat[,id]
  
  # filter modules strongly correlated with Diagnosis and cell type markers
  cellres = read.csv("./results/tables/network/CellTypeEnrichmentFDR.csv")
  cellres = melt(cellres)
  cellres$variable = gsub(".FDR","",cellres$variable)
  cellres = subset(cellres, value<0.05)
  
  ## load correlation with covariates 
  load("./data/network/mod_cov_cor.Rdata")
  cordat = melt(moduleCors)
  colnames(cordat) = c("Dx","mod","cor")
  
  p = melt(modulesPval)
  cordat$p = p$value
  cordat = subset(cordat,Dx == "Diagnosis" & p<0.05)
  idx_mod = union(cordat$mod,cellres$variable)
  
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
  df = df[df$m %in% idx_mod,]
  
  df$symbol = ""
  df$symbol[df$p<0.05] = "*"
  
  # df$cell = celldat$cell[match(df$module,celldat$module)]
  
  col = df$color
  names(col) = df$m
  df$module = factor(df$module,levels = unique(df$module))
  df$diagnosis = factor(df$diagnosis,levels = c("AsymAD","AD"))
  
  ## filter modues with signifiant anova
  anov_idx = unique(df$module[df$anova<0.05])
  pdat = df[df$module %in% anov_idx ,]
  txt = pdat[pdat$anova<0.05 & pdat$diagnosis =="AD",]
  txt$anova = signif(txt$anova,2)
  
  # plot
  p = ggplot(pdat, aes(x = diagnosis,
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
              aes(y = beta + sign(beta) * se + sign(beta) * .001),
              position = position_dodge(width = 0.5)) +
    geom_text(data = txt,aes(label = paste0("P = ",anova),
                             y = 0.075, x= 1.5),size=3)+
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
         filename = "./results/figures/network/sex_cell_mod_beta.pdf",
         width = 14,
         height = 7)
  
  ## correlation for sex-specific modules
  idx = c("M4","M10","M16","M20","M25")
  dat = data.frame()
  alldat = cbind(eigmat[,idx],datMeta)
  res = data.frame()
  
  for(s in sex){
    dat = alldat[alldat$Sex ==s,]
    
    for(m in idx){
      r = cor.test(dat[,m],dat[,"Braak"])
      rho = unname(r$estimate)
      p = unname(r$p.value)
      res = rbind(res,cbind(rho,p,m,s))
    }
  }
  
  res$rho = as.numeric(res$rho)
  res$p = as.numeric(res$p)
  df = tidyr::pivot_longer(alldat,cols = 1:5,names_to = "module",values_to = "eigenprotein")
  
  ggplot(df,aes(x = ABETA.2pep.corrected,y = eigenprotein,color= Sex))+
    geom_point()+
    facet_wrap(~Sex)+
    geom_smooth(color = "black",se = T)
  
}
