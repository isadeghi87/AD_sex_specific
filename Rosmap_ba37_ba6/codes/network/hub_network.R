## hub_network.R
### gene network of top hub genes for M7
rm(list=ls()); options(stringsAsFactors = F)
pacman::p_load(gplots,WGCNA,ggpubr,dplyr,magrittr,
               tidyverse,ggthemes,igraph)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
load("./data/network/parameters/finalNetwork.RData")
## hub genes
kME = signedKME(t(datExp), MEs$eigengenes)
colnames(kME) = gsub("kME", "M", colnames(kME))

n = 50  #plot top 5 hub genes for each module
res = data.frame()
for(i in colnames(kME)){
  genes_idx = order(kME[which(mods == i),i], decreasing = T)[1:n]
  hubs = rownames(datExp[genes_idx,])
  # symb = gene[match(hubs,rownames(datExp))]
  res = rbind(res, cbind(i,hubs))
}

res = res[res$i=="M7",]
res$symb = tableSup$gene[match(res$hubs,rownames(tableSup))]
keepgenes = rownames(tableSup)[tableSup$gene %in% res$symb]
gene_idx = res$hubs

# keep hubgenes
adjMat = bicor(t(datExp))
adjMat = adjMat[gene_idx,gene_idx]

## keep edges with corr > 0.65
topcors=0.85^16
adjMat[adjMat< topcors] = 0

g1 <- graph.adjacency(as.matrix(adjMat),
                      mode="undirected",
                      weighted=T,
                      diag=FALSE)

mds = cmdscale(dist(t(adjMat)), eig = T)
layoutFR = mds$points
mdsedgecolors = numbers2colors(E(g1)$weight, 
                               colors = rev(blueWhiteRed(100, gamma=1)),
                               signed=T,
                               centered=T,
                               lim=c(-1,1))

# remove weigth edges < 0.7
g.copy <- delete.edges(g1, which(E(g1)$weight <0.2))

pdf("./results/figures/network/TopHubGenes.pdf", width=5, height=4)
par(mar=c(1,1,1,1))
set.seed(123)
plot.igraph(g.copy, 
            vertex.label = res$symb,
            vertex.label.dist= 1, 
            edge.width=E(g1)$weight,
            vertex.size=5, 
            vertex.frame.color="black",
            vertex.label.color="black",
            vertex.color = "black",
            vertex.label.cex=0.4, 
            layout=layout.fruchterman.reingold(g1),
            edge.color=mdsedgecolors)
title(main = "", cex.main = 1)

dev.off()
