#st5_network_module_DE.R
if(T){
rm(list=ls()); options(stringsAsFactors = F)
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme);library(dplyr)
library(magrittr);library(ggthemes)
condition = TRUE

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
sex = c("female","male")

cell_modules = readRDS("./codes/network/cell_modules.rds")
plot.list = list()  
res = data.frame()

for (s in sex){
  load(paste0("./codes/network/finalNetwork_",s,".RData",sep=""))
  celldat = subset(cell_modules,sex ==s)
  #---load final network
  eigmat = MEs$eigengenes
  colnames(eigmat) = gsub("ME","M",colnames(eigmat))
  eigmat = eigmat[,-1]
  kME = signedKME(t(exp), MEs$eigengenes)
  colnames(kME) = gsub("kME", "M", colnames(kME))
  kME = kME[,-1]
  
all_colors = colnames(eigmat)

#Step 1 - Calculate module-trait P values, beta, and SEM
moduleTraitP = matrix(NA, nrow = length(all_colors),
                      ncol=3)
colnames(moduleTraitP) = c("ANOVA","AD","PA")
rownames(moduleTraitP) = all_colors
moduleTraitB = moduleTraitSE = moduleTraitP

if(T){
  for (m in all_colors) {
    me = eigmat[[m]]
    
    model = lm(me ~ Dx + RIN , data = meta)
    moduleTraitP[m,"ANOVA"] = anova(model)["Dx","Pr(>F)"]
    mixedmodel = summary(model)
    mixedmodel = mixedmodel$coefficients
    for(var in c("AD", "PA")) {
      moduleTraitB[m,var] = mixedmodel[grep(var, rownames(mixedmodel)), 1]
      moduleTraitSE[m,var] = mixedmodel[grep(var, rownames(mixedmodel)), 2]
      moduleTraitP[m,var] = mixedmodel[grep(var, rownames(mixedmodel)), 4]
    }
  }
}
  moduleTraitP.fdr = p.adjust(moduleTraitP, "fdr")
  dim(moduleTraitP.fdr) = dim(moduleTraitP)
  dimnames(moduleTraitP.fdr) = dimnames(moduleTraitP)
  
  ## plot for all modules
  bpdata = melt(moduleTraitB[,-1])
  sedata = melt(moduleTraitSE[,-1])
  pdata = melt(moduleTraitP[,-1])
  bpdata$se = sedata$value
  bpdata$p = pdata$value
  
  colnames(bpdata)[1:3] = c("module","condition","beta")
  bpdata$symbol = ""
  bpdata$symbol[bpdata$p<0.05] = "*"
  bpdata$symbol[bpdata$p<0.01] = "**"
  bpdata$symbol[bpdata$p<0.001] = "***"
  bpdata$cell = celldat$cell[match(bpdata$module,celldat$module)]
  bpdata$color = celldat$color[match(bpdata$module,celldat$module)]
  bpdata$sex = s
  res = rbind(res,bpdata)
  

res = na.omit(res)
cols = unique(res$color)
names(cols) = res$module[match(cols,res$color)]

# plot
p = ggplot(res, aes(x = condition,
                       y = beta,
                      group = module,
                       fill = module,
                       label = symbol)) + 
  facet_grid( sex~ cell) +
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
  scale_fill_manual(values = cols) +
  labs(y = "beta", x = "", title = "") +
  scale_x_discrete()+
  theme_few(base_size = 14)+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 14,face = "bold",color = "white"),
        axis.title = element_text(size = 12,angle = 90),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background.y = element_rect(fill = "darkblue"),
        strip.background.x = element_rect(fill = "darkred"),
        )

p
ggsave(plot = p,
       filename = "./results/figures/network/sex_cell_mod_beta.pdf",
       width = 10,
       height = 5)
}
}
