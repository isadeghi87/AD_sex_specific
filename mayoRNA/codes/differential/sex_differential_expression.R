### sex_differential_expression.R
## Here we calculate sex-specific differential expression 
if(T){
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/")
  
  pacman::p_load(reshape,limma,ggplot2,readr,dplyr,purrr,ggpubr,ggthemes,
                 edgeR,stringr,tidyr,statmod)
  
  
  #### load data from bulk analysis ####
  load("./codes/covariates/normalized_data.Rdata")
  
  #### differential expression ####
  sex = unique(datmeta$Sex)
  
  stats = data.frame()
  # stratify for sex
  for (s in sex) {
    idx = which(datmeta$Sex == s)
    Meta = datmeta[idx,]
    Expr = datExp[,idx]
    mod <- model.matrix(~ Dx+ Age + RIN, data = Meta)
    colnames(mod) = gsub("Dx","",colnames(mod))
    fit <- lmFit(Expr, mod)
    
    efit <- eBayes(fit = fit, trend = T, robust = T)
    ## PA stats
    meta_pa <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
    meta_pa <- meta_pa[,c(1,4,5)]
    meta_pa$sex = s
    meta_pa$Dx = "PA"
    meta_pa$gene_id = rownames(meta_pa)
    rownames(meta_pa) = NULL
    
    ## AD stats
    meta_ad <- topTable(efit, coef = 3, number = Inf, sort.by = "none")
    meta_ad <- meta_ad[,c(1,4,5)]
    meta_ad$sex = s
    meta_ad$Dx = "AD"
    meta_ad$gene_id = rownames(meta_ad)
    rownames(meta_ad) = NULL
    
    stats = rbind(stats,meta_pa,meta_ad)
    
  }
  
  colnames(stats)[2:3] = c("p.value","fdr") 
  write.csv(stats, "./results/tables/AD_PA_sex_sumstats.csv")
  
  # spread the data
  sig = subset(stats, fdr < 0.05 & abs(logFC)>0.58)
  table(sig$Dx,sig$sex)
  
  ### volcano plot 
  p = ggplot(stats, aes(x = logFC,y = -log10(fdr)))+
    geom_point(fill = "grey",color = "black", shape = 21, size = 2)+
    geom_point(data=sig, aes(x=logFC, y = -log10(fdr)),fill = "blue",
               shape = 21, size = 2,color = "black")+
    # scale_color_manual(values = c("green","red"))+
    facet_grid(Dx ~ sex)+
    labs(title = "Differentially expressed genes",
         subtitle = "FDR <0.05")+
    theme_bw(base_size = 14)+
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "black"),
          strip.text = element_text(face="bold",size=14,color= "white"))
  
  p
  
  ggsave("./results/figures/differential/volcano_plot.pdf",plot=p,
         width = 8,height = 5.5)
  
  save(file = "./codes/differential/sex_differential_expression.Rdata",sig,stats)
}
