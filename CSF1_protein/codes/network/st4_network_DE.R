#st5_network_module_DE.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
 library(ggplot2); 
  library(reshape); library(nlme);library(dplyr)
  library(magrittr);library(ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF1/")
  sex = c("female","male")
  
  #load matched modules ####
  matched = read.csv("../results/tables/network/preservation/matched_modules_results.csv")
  matched = matched[which(matched$rep !="M0"),]
  matched = matched[matched$dataset =="CSF1",]
  res = data.frame()
  
  for (s in sex){
    celldat = subset(matched,sex ==s)
    id_match = celldat$rep
    
    #---load final network
    load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
    eigmat = MEs$eigengenes
    colnames(eigmat) = gsub("ME","M",colnames(eigmat))
    id = which(colnames(eigmat)=="M0")
    eigmat = eigmat[,-id]
    
    #Step 1 - Calculate module-trait P values, logFC, and SEM
    df = as.data.frame(matrix(NA, nrow = ncol(eigmat),
                                ncol=3))
    colnames(df) = c("beta","se","p")
    rownames(df) = colnames(eigmat)
    
    for (m in colnames(eigmat)) {
      me = as.numeric(eigmat[[m]])
      
      mod = lm(me ~ Diagnosis, data = meta)
      mod = summary(mod)$coefficients
      df[m,1] = mod[2, 1]# beta
      df[m,2] = mod[2, 2]# std.error
      df[m,3] = mod[2, 4]# p value
    }
    
    df$fdr = p.adjust(df$p, "fdr")
    df$color = names(mods)[match(rownames(df),mods)]
    df$module = rownames(df)
    df = df[df$module %in% id_match,]
    df$symbol = ""
    df$symbol[df$p<0.05] = "*"
    df$symbol[df$p<0.01] = "**"
    df$symbol[df$p<0.001] = "***"
    df$symbol[df$p<0.001] = "***"
    
    # df$cell = celldat$cell[match(df$module,celldat$module)]
    df$sex = s
    res = rbind(res,df)
  }
  
  col = res$color
  names(col) = res$module
  res$module = factor(res$module,levels = unique(res$module))
  # plot
  p = ggplot(res, aes(x = module,
                      y = beta,
                      fill = module,
                      label = symbol)) + 
    facet_wrap(~sex, scales = "free",ncol =1) +
    geom_bar(stat = "identity",
             width = 0.5,
             position = position_dodge(),
             color = "black") +
    geom_errorbar( aes(ymin = (beta - se), ymax = (beta + se)),
                   position = position_dodge(width = 0.5),
                   # size = 0.2,
                   width = 0.2) +
    geom_text(color = "red",
              size = 4,
              aes(y = beta + sign(beta) * se + sign(beta) * .0009),
              position = position_dodge(width = 0.5)) +
    scale_fill_manual(values = col) +
    labs(y = "beta", x = "", title = "CSF1 modules expression") +
    scale_x_discrete()+
    theme_bw(base_size = 14)+
    theme(legend.position = "none",
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 14,face = "bold",color = "white"),
          axis.title = element_text(size = 12,angle = 90),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.background.y = element_rect(fill = "darkblue"),
          strip.background.x = element_rect(fill = "darkgreen"))
  
  p
    ggsave(plot = p,
           filename = "./results/figures/network/sex_cell_mod_beta.pdf",
           width = 5,
           height = 5)
}
