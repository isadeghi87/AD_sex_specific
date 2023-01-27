#5c_network_moduleTrait.R
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme);library(dplyr)
library(magrittr);library(ggthemes)
condition = TRUE

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/bulk_analysis/replication_msbb/")

# load final network
load("./codes/network/parameters/finalNetwork.RData")

eigmat = MEs$eigengenes
eigmat= eigmat[,-1]
colnames(eigmat) = gsub("E","",colnames(eigmat))

all_colors = colnames(eigmat)
all_genes = colors
names(all_genes) = rownames(datExp)

#Step 1 - Calculate module-trait P values, beta, and SEM
moduleTrait = as.data.frame(matrix(NA, nrow = length(all_colors),
                      ncol=3))
colnames(moduleTrait) = c("beta","SE","p")
rownames(moduleTrait) = all_colors

if(condition){
  for (m in all_colors) {
    me = eigmat[[m]]
    
    mixedmodel = lme(me ~ Dx , data = datMeta, random = ~1|Subject_ID)
    mixedmodel = summary(mixedmodel)$tTable
      moduleTrait[m,1] = mixedmodel[2, 1]# beta
      moduleTrait[m,2] = mixedmodel[2, 2]# std.error
      moduleTrait[m,3] = mixedmodel[2, 5]# p value
    }
}

moduleTrait$fdr = p.adjust(moduleTrait$p, "fdr")
moduleTrait$color = names(mods)[match(rownames(moduleTrait),mods)]
moduleTrait$module = rownames(moduleTrait)

moduleTrait$symbol = ""
moduleTrait$symbol[moduleTrait$p<0.05] = "*"
moduleTrait$symbol[moduleTrait$p<0.01] = "**"
moduleTrait$symbol[moduleTrait$p<0.001] = "***"

## assign colors
mod_col = moduleTrait$color
names(mod_col) = moduleTrait$module
  
p = ggplot(moduleTrait, aes(x = module,
                         y = beta,
                         fill = module,
                         label = symbol)) + 
    geom_bar(stat = "identity",width = 0.7,
             position = position_dodge(),
             color = "black") +
    geom_errorbar(aes(ymin = (beta - SE), ymax = (beta + SE)),
                   position = position_dodge(.9),
                   size = 0.3,
                   width = .3) +
    geom_text(color = "red",
              size = 4,
              aes(y = beta + sign(beta) * SE + sign(beta) * .0005),
              position = position_dodge(0.1)) +
    scale_fill_manual(values = mod_col) +
    labs(y = "beta", x = "", title = expression("Transcriptome replication"[msbb]*" modules expression")) +
    scale_x_discrete()+
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 12,angle = 90),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "#0A043C"))

ggsave("./results/figures/network/msbb_module_DE_supp.pdf",width = 6,height = 3,plot = p)

### keep only cell type modules
neuron = c("M1","M6")
endo = "M3"
astro = "M4"
oligo = "M5"
micro = "M10"
cellmod = c(neuron,astro,oligo,micro,endo)

moduleTrait = moduleTrait[cellmod,]
moduleTrait$cell = moduleTrait$module
moduleTrait$cell[moduleTrait$module == micro] = "microglia"
moduleTrait$cell[moduleTrait$module %in% neuron] = "neuron"
moduleTrait$cell[moduleTrait$module == astro] = "astrocyte"
moduleTrait$cell[moduleTrait$module == oligo] = "oligo"
moduleTrait$cell[moduleTrait$module == endo] = "endothelial"
moduleTrait$cell = factor(moduleTrait$cell,
                          levels =unique(moduleTrait$cell))

## order of colors for modules
mod_col = moduleTrait$color
names(mod_col) = moduleTrait$module 

p = ggplot(moduleTrait, aes(x = module,
                            y = beta,
                            fill = module,
                            label = symbol)) + 
  geom_bar(stat = "identity",width = 0.3,
           position = position_dodge(),
           color = "black") +
  geom_errorbar(aes(ymin = (beta - SE), ymax = (beta + SE)),
                position = position_dodge(.9),
                size = 0.2,
                width = .2) +
  geom_text(color = "red",
            size = 4,
            aes(y = beta + sign(beta) * SE + sign(beta) * .0001),
            position = position_dodge(0.9)) +
  scale_fill_manual(values = mod_col) +
  facet_wrap(~cell,nrow =1 ,scales = "free_x")+
  labs(y = "beta", x = "", title = expression("Replication"[msbb]*" modules expression")) +
  scale_x_discrete()+
  theme_bw()+
  theme(
    legend.position = "none",
    strip.text.x = element_text(color = "white", face = "bold"),
    axis.title.y = element_text(size = 12,angle = 90),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold",size = 14),
    strip.background = element_rect(fill = "#0A043C"))
p

ggsave(plot = p,
       filename = "./results/figures/network/msbb_module_DE.pdf",
       width = 7,
       height = 3)


save.image("./codes/network/mod_DE.Rdata")
