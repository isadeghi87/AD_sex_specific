#st5_networkVisualization.R

rm(list=ls())
options(stringsAsFactors = F)
#source("http://bioconductor.org/biocLite.R"); biocLite("igraph")
library(WGCNA);library(ggplot2); library(reshape); 
library(igraph); library(RColorBrewer); 
library(WGCNA); library(corrplot);library(ggthemes)
library(dplyr);library(corrr);library(igraph)

condition = TRUE

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/")

#load final network
load("./codes/wgcna/finalNetwork.RData")

#### Make module eigengene-MDS plot ####
eigmat = MEs$eigengenes
colnames(eigmat) = gsub("E","",colnames(eigmat))
eigmat = eigmat[,-1]
adj = bicor(eigmat)
cols = cols[-1]

#### correlation between modules ####

pdf(file = "./results/wgcna/modules_cor.pdf",
    height = 10, 
    width = 10)
corrplot(signif(adj,1),
         # col = blueWhiteRed(n = 1000, gamma = 1)[1000:1],
         # main = "Gene modules pairwise correlation",
         mar = c(0,0,0,0),
         method = "ellipse",
         type = "upper",
         addCoef.col = T,
         order = "hclust",
         hclust.method = "complete",
         tl.cex = 1.5,
         tl.srt = 45,
         is.corr = T,
         diag = F,
         tl.col = "black",
         number.cex = 1.2)
dev.off()

net = network_plot(adj,curved = F,min_cor = 0.5)
ggsave("./results/wgcna/mod_network_cor.pdf",
       width = 8,height = 6, plot = net)



## Part 2) hub genes #### 
mod_gene = data.frame(module = mods,
                       row.names = rownames(datExp))

cons_kme = signedKME(t(datExp), MEs$eigengenes)
cons_kme = cons_kme[,-which(colnames(cons_kme)=="kME0")]
colnames(cons_kme) <- gsub("kME", "M", colnames(cons_kme))

maxsize = 10  #plot top 20 hub genes for each module
gene_idx = order(cons_kme[which(colors=="turquoise"),1], decreasing = T)[1:maxsize]

for(i in c(2:ncol(cons_kme))){
  gene_idx = c(gene_idx, order(cons_kme[,i], decreasing = T)[1:maxsize])
}

hubGenes = character()
hubGenes.kme = numeric()

#---hub for cell-type modules
oligo = 'M1'
microglia = 'M4'
neuron  = c('M5','M9','M11','M15')
astro = c('M13','M15')
endo = 'M6'
cellMod = c(oligo,microglia,neuron,astro,endo)

for(col in cellMod)  {
  modgenes = rownames(datExp)[which((mods) == col)]
  kmes = cons_kme[modgenes, col]
  top_hubs = modgenes[order(kmes, decreasing=T)[1:maxsize]]
  top_hubs.kme = kmes[order(kmes, decreasing=T)[1:maxsize]]
  hubGenes = c(hubGenes,top_hubs)
  hubGenes.kme = c(hubGenes.kme, top_hubs.kme)
}
gene_idx = match(hubGenes,rownames(datExp))

#### save hubgenes
saveRDS(object = hubGenes,"./tables/wgcna/hubGenes.rds")

## calculate adjacency matrix 
adjMat = bicor(t(datExp))
keepgenes = rownames(cons_kme)[gene_idx]

# keep hubgenes
adjMat = adjMat[gene_idx,gene_idx]

## keep edges with corr > 0.65
topcors=0.3
adjMat[adjMat< topcors] = 0

geneSymbols = datProbes$external_gene_name[match(keepgenes, datProbes$ensembl_gene_id)]
g1 <- graph.adjacency(as.matrix(adjMat),
                      mode="undirected",
                      weighted=T,
                      diag=FALSE)

mds = cmdscale(dist(t(adjMat)), eig = T)
layoutFR = mds$points
mdsedgecolors = numbers2colors(E(g1)$weight, 
                               colors = blueWhiteRed(100, gamma=1),
                               signed=T,
                               centered=T,
                               lim=c(-1,1))

# remove weigth edges < 0.7
g.copy <- delete.edges(g1, which(E(g1)$weight <0.9))
mst <- mst(g1, algorithm="Prim")
# mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
# mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
# V(mst)$color <- mst.communities$membership+1 

pdf("./results/wgcna/TopHubGenes_cells.pdf", width=5, height=5)
par(mar=c(0,0,0,0))
set.seed(123)
plot(g.copy, 
     vertex.label = geneSymbols,
     vertex.label.dist=0.7, 
            edge.width=0.1,
            vertex.size=4, 
            vertex.frame.color="black",
            vertex.label.color="black",
            vertex.color = colors[gene_idx],
            vertex.label.cex=0.4, 
            layout=layout.fruchterman.reingold(g1),
            edge.color="grey")
title(main = "", cex.main = 1)

dev.off()

save.image("./codes/wgcna/St5_networkVis.Rdata")
