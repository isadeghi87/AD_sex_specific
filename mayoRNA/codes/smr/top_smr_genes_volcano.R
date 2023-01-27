## expression of top genes from mendelian randomization (SMR)
if(T){
rm(list=ls())
library(dplyr); 
library(reshape2);library(ggplot2);
library(ggrepel);library(ggthemes);library(ggplot2)

## sumstats of DGE
stat = read.csv("./results/tables/AD_PA_sex_sumstats.csv") 
## results of SMR
topsnp = read.csv("./results/tables/network/topSNP_smr.csv")
rs = stat[stat$gene_id %in% unique(topsnp$probeID),]
rs$gene = topsnp$Gene[match(rs$gene_id,topsnp$probeID)]

p = ggplot(rs, aes(x = logFC, y = -log10(fdr)))+
  geom_point(aes(fill = sex),color = "black",shape=21,size=3,show.legend = F)+
  facet_grid(Dx~sex)+
  # scale_fill_manual(values = c("magenta","blue"))+
  geom_hline(yintercept = -log10(0.05),lty=2,color = "red")+
  geom_text_repel(data = rs[rs$fdr<0.05,],aes(label = gene),
                  size = 2,color = "blue")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "darkblue"),
        strip.text = element_text(colour = "white",face = "bold",size=14))

ggsave(plot=p, "./results/figures/smr/top_SMR_genes_volcano.pdf",
       width = 6,height = 4) 
}