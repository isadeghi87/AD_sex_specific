#5c_network_meta.R
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme);library(dplyr)
library(magrittr);library(ggthemes)
condition = TRUE

SEtwd("C:/USErs/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/diSEaSE_specific/AD_PA_CTL/bulk_analysis/replication_proteomics_CSF1/")

# load final network
load("./codes/network/parameters/finalNetwork.RData")

eigmat = MEs$eigengenes
eigmat= eigmat[,-1]
colnames(eigmat) = gsub("E","",colnames(eigmat))
all_colors = colnames(eigmat)

meta = as.data.frame(matrix(NA, nrow = length(all_colors),
                            ncol=3))
colnames(meta) = c("logFC","SE","p")
rownames(meta) = all_colors

if(condition){
  for (m in all_colors) {
    me = eigmat[[m]]
    
    mod = lm(me ~ Diagnosis  , data = datMeta)
    mod = summary(mod)$coefficients
    meta[m,1] = mod[2, 1]# logFC
    meta[m,2] = mod[2, 2]# std.error
    meta[m,3] = mod[2, 4]# p value
  }
}

meta$fdr = p.adjust(meta$p, "fdr")
meta$color = names(mods)[match(rownames(meta),mods)]
meta$module = rownames(meta)

meta$symbol = ""
meta$symbol[meta$fdr<0.05] = "*"
meta$symbol[meta$fdr<0.01] = "**"
meta$symbol[meta$fdr<0.001] = "***"
meta$symbol[meta$fdr<0.001] = "***"
meta$module = factor(meta$module,levels = unique(meta$module))

## assign colors
mod_col = meta$color
names(mod_col) = meta$module

p = ggplot(meta, aes(x = module,
                       y = logFC,
                       fill = module,
                       label = symbol)) + 
  geom_bar(stat = "identity",width = 0.3,
           position = position_dodge(),
           color = "black") +
  geom_errorbar(aes(ymin = (logFC - SE), ymax = (logFC + SE)),
                position = position_dodge(.9),
                size = 0.2,
                width = .2) +
  geom_text(color = "red",
            size = 4,
            aes(y = logFC + sign(logFC) * SE + sign(logFC) * .0001),
            position = position_dodge(0.9)) +
  scale_fill_manual(values = mod_col) +
  # facet_wrap(~module,ncol = 3,scales = "free")+
  labs(y = "logFC", x = "", 
       title = expression("Replication proteomics"[CSF]*" modules expression")) +
  scale_x_discrete()+
  theme_bw()+
  theme(
    legend.position = "none",
    strip.text.x = element_text(color = "white",face = "bold"),
    axis.title.y = element_text(size = 12,angle = 90),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold",size = 14),
    strip.background = element_rect(fill = "#0A043C"))

p
ggsave(p,filename = "./results/figures/network/csfProteomics_modules_DE_supp.pdf",
       width = 8,height = 6)

### keep only cell type modules
neuron = "M1"
astro.neuro = "M6"
astro.olig = "M2"
micro = c("M3","M4")
cellmod = c(neuron,astro.neuro,astro.olig,micro)

meta = meta[meta$module %in% cellmod, ]
meta$cell = as.character(meta$module)
meta$cell[meta$module %in% neuron] = "neuron"
meta$cell[meta$module == astro.micro] = "astro/micro"
meta$cell[meta$module == micro] = "microglia"
meta$cell = factor(meta$cell,
                   levels = unique(meta$cell))


p = ggplot(meta, aes(x = module,
                       y = logFC,
                       fill = module,
                       label = symbol)) + 
  geom_bar(stat = "identity",width = 0.5,
           position = position_dodge(0.9),
           color = "black") +
  geom_errorbar(aes(ymin = (logFC - SE), ymax = (logFC + SE)),
                position = position_dodge(.9),
                size = 0.1,
                width = .2) +
  geom_text(color = "red",
            size = 4,
            aes(y = logFC + sign(logFC) * SE + sign(logFC) * .0005),
            position = position_dodge(0.9)) +
  scale_fill_manual(values = mod_col) +
  # facet_wrap(~ cell,nrow = 1,scales = "free_x")+
  labs(y = "logFC", x = "", 
       title = expression("Replication proteomics"[CSF]*" modules expression")) +
  # scale_x_discrete()+
  theme_bw()+
  theme(
    legend.position = "none",
    strip.text.x = element_text(color = "white",face = "bold"),
    axis.title.y = element_text(size = 12,angle = 90),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold",size = 14),
    strip.background = element_rect(fill = "#0A043C"))

p
ggsave(plot = p,
       filename = "./results/figures/network/csfProteomics_cell_modules_DE.pdf",
       width = 6,
       height = 4)


save.image("./codes/network/mod_DE.Rdata")
