## st8_AD_genes_enrichment.R
if(T){
rm(list=ls()); options(stringsAsFactors = F)
library(gplots); library(biomaRt); library(WGCNA)
library(ggpubr); library(dplyr); library(magrittr)
library(GeneOverlap)
library(readxl);library(tidyverse);library(ComplexHeatmap)
condition = TRUE

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
source("../bulk_analysis/codes/source/Fisher_overlap.R")

sex = c("female","male")
# cell-type modules: 
cellmod = readRDS("./codes/network/cell_modules.rds")
enrich = data.frame()

### list of AD risk factor genes
ad_dat = readxl::read_xlsx("../data/gwas/gwas.xlsx",
                           sheet = 20, skip =5) 
adGene= unique(ad_dat$`GWAS Gene`)

sex = c("female","male")
for(s in sex){
  #load final network
  load(paste0("./codes/network/finalNetwork_",s,".RData",sep=""))
  
eigmat = MEs$eigengenes
eigmat = eigmat[,-1]
colnames(eigmat) = gsub("ME","M",colnames(eigmat))
kME = signedKME(t(exp), MEs$eigengenes[,-1])
colnames(kME) = gsub("kME", "M", colnames(kME))

## check enrichment ####
res = data.frame(module = "NA",enrichment="NA",p="NA")

for (m in colnames(kME)) {
  id1 = datProbes$gene_name[mods==m]
  f = newGeneOverlap(id1,
            adGene,
            spec = "hg19.gene" )
    obj = testGeneOverlap(f)
    en = signif(getOddsRatio(obj),2)
    p = signif(getPval(obj),2)
    res = rbind(res,c(m,en,p))
}

res = res[-1,]
res$FDR = p.adjust(res$p,"fdr")
res$color = names(mods)[match(res$module,mods)]
res$module = factor(res$module,levels = res$module)
res$symbol = ""
res$symbol[res$FDR<0.05]="*"
res$symbol[res$FDR<0.01]="**"
res$symbol[res$FDR<0.001]="***"
res$enrichment = as.numeric(res$enrichment)
res$FDR = as.numeric(res$FDR)
res$sex = s
enrich = rbind(enrich,res)

}
}

## draw plot
col = ifelse(res$FDR<0.05,"blue","white")
colManuel = names(mods)
names(colManuel) = mods

p   = ggplot(enrich, aes(x = module,
                      y = enrichment,
                      # alpha = FDR<0.05,
             fill = module,
             label = symbol))+
  geom_bar(stat = "identity",
             colour = "black",width = 0.5 )+
  facet_wrap(~sex,ncol=1,scales = "free")+
  geom_text(aes(label = symbol, y = enrichment +0.1), color = "red",size = 5, show.legend = F)+
  scale_fill_manual(values = colManuel)+
  labs(title = "Enrichment of late-onset AD genes",
       alpha = "Significant\n(FDR<0.05)", x = "")+
  theme_bw(base_size = 14)+
  theme(axis.text.x.bottom = element_text(angle= 45,vjust = 0.5),
        plot.title = element_text(hjust=0.5, size =  14),
        strip.background = element_rect(fill = "darkblue"),
        strip.text = element_text(size = 14,face = "bold",color = "white"))+
  guides(fill = "none", size = "none")
p
ggsave("./results/figures/network/AD_gene_enrichment.pdf",width = 10,height = 5,plot = p)

