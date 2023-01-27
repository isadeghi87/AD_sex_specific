## we look at the expression of top important genes from
# classifier models
if(T){
  rm(list=ls())
  library(dplyr); 
  library(reshape2);library(ggplot2);
  library(ggrepel);library(ggthemes);library(ggplot2)
 
   setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  sex = c("female","male")
  
  stat = read.csv("./results/tables/AD_PA_sex_sumstats.csv")
  gene = read.csv("./results/tables/conditions_top_important_genes.csv")
  # id = unique(gene$var)
  # df = stat[stat$gene_id %in% id,]
  # df$gene_name = gene$gene_name[match(df$gene_id,gene$var)]
  # 
  res = data.frame()
  for(s in sex){
    dat = stat[stat$sex==s,]
    id = gene$var[gene$s==s]
    dat = dat[dat$gene_id %in%  id,]
    dat$gene_name = gene$gene_name[match(dat$gene_id,gene$var)]
    
    res = rbind(res,dat)
  }
  
  p =  ggplot(res, aes(x = logFC, y = -log10(fdr)))+
    geom_point(aes(fill = sex),color = "black",
               shape=21,size=2,show.legend = F)+
    facet_grid(Dx~sex)+
    geom_hline(yintercept = -log10(0.05),lty=2,color = "red")+
    geom_text_repel(data = res[res$fdr<0.05,],aes(label = gene_name),
                    size = 2,color = "blue")+
    theme_bw()+
    theme(strip.background = element_rect(fill = "darkblue"),
          strip.text = element_text(colour = "white",face = "bold",size=14))
  
  ggsave(plot=p, "./results/figures/models/models_important_genes_volcano.pdf",
         width = 7,height = 5) 
}