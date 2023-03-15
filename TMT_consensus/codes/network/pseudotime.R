setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
load("./data/network/parameters/finalNetwork.RData")
source("C:/Users/crgcomu/Desktop/Iman/Brain_meta/source/lineageFunctions.R")
library(monocle)

#Performing lineage inference with Monocle2
id = which(tableSup$mods=="M7")
temp <- datExp[id,]
sex = c("female","male")
ps.plot = dx.plot = brak.plot = list()

for(s in sex){
id = which(datMeta$Sex==s)
dat = temp[,id]
meta = datMeta[id,]
rownames(dat) <- NULL
colnames(dat) <- NULL
MonRun <- RunMonocleTobit(dat, meta$Diagnosis)

meta$pseudotime = MonRun$Pseudotime
meta$Braak =factor(meta$Braak)
ps = ggplot(meta,aes(x = Diagnosis,y = pseudotime,fill = Diagnosis))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Dark2")+
  labs(x = "",title=s)
ps.plot[[s]] = ps

## Visualizing using Monocle's visualization 
# dx = plot_cell_trajectory(MonRun, color_by = meta$Diagnosis)
# brak = plot_cell_trajectory(MonRun, color_by = meta$Braak,show_backbone = T,show_tree = T,
#                           show_cell_names = F,show_branch_points = F)

## Visualizing Monocle2 with different labels 

dx = Monocle.Plot(MonRun, Labels = meta$Diagnosis, Discrete = T)+
  labs(title = s,x = "component 1",y = "component 2")+
  scale_color_brewer(palette = "Dark2")+
  geom_smooth(color = "black")+
  guides(fill = guide_legend(title="Diagnosis"))+
 theme_bw()

brak = Monocle.Plot(MonRun, Labels = meta$Braak, Discrete = F)+
  labs(title = s,x = "component 1",y = "component 2")+
  scale_color_continuous(type = "viridis")+
  guides(color = guide_legend(title="Braak stage"))+
  theme_bw()


dx.plot[[s]] = dx
brak.plot[[s]] = brak
}

library(ggpubr)
p1 = ggarrange(plotlist = dx.plot,legend = "right",common.legend = T)
p2 = ggarrange(plotlist = brak.plot,legend = "right",common.legend = T)

gg = ggarrange(p1,p2,ncol = 1)
ggsave(filename = "./results/figures/network/pseudotime.pdf",gg,
       width = 8,height = 6)
psplot = ggarrange(plotlist = ps.plot,common.legend = T)

ggsave(filename = "./results/figures/network/pseudoBoxPlot.pdf",psplot,
       width = 5,height = 3)

