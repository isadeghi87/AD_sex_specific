### overlap.R
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/")

## libraries

library(ggplot2);library(RRHO);library(tidyverse);
library(corrr);library(corrplot);library(GGally)
library(network);library(sna);library(ggplot2)

### load data
load("./codes/differential_expression.Rdata")
sumstats = read.csv("./tables/AD_PA_sumstats.csv",row.names = 1)
sumstats = select(sumstats,contains("logFC"))
colnames(sumstats) = gsub(".logFC","",colnames(sumstats))

### correlation
library(ggstatsplot)
library(ggcorrplot)
set.seed(123)

# plot
cor_p = ggcorrmat(
  data = sumstats,p.adjust.method = "fdr",
  ggtheme = theme_bw(base_size = 10),
  type = "spearman", # correlation method
  colors = c("#FF449F", "white", "#005F99"),
  matrix.type = "lower" # type of matrix
)

ggsave("./results/overlap/corr.pdf",width = 8,height = 6,cor_p)

#####
# random graph
cor_dat = cor(sumstats,use = "pairwise.complete.obs",method = "spearman")

### corr network
library(corrr)

p = network_plot(cor_data,min_cor = 0,curve=T,legend= T)
ggsave(filename = "./results/overlap/Correlation_network.pdf",
       plot = p,width = 8,height = 4)

#### RRHO plots ####

rrho.list <- list()
for ( i in 1:ncol(sumstats)){
  df <- data.frame(gene = rownames(sumstats),
                   logfc = sumstats[,i])
  rrho.list[[i]] = df
}

names(rrho.list) <- colnames(sumstats)

jet.colors  <- colorRampPalette(
  c("#00007F", "blue", "#007FFF", "cyan", 
    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));

plot = T
if(plot){
pdf("./results/overlap/rrho_sex_comparisons.pdf", width = 3, height = 12)
par(mfrow=c(4,1))
# par(mar = c(0.2, 0.2, 0, 0))

for (i in c(1,3,5,7)){
  R.obj <- RRHO(list1 = rrho.list[[i]], 
                list2 = rrho.list[[i+1]],
                BY = T,
                alternative = "enrichment",
                plots = T,
                # log10.ind = T,
                outputdir = "./results/overlap/rrho/", 
                labels = c(names(rrho.list)[i], names(rrho.list)[i+1])
  )
  
  image(R.obj$hypermat,
        xlab = names(rrho.list)[i],
        ylab=names(rrho.list)[i+1],
        col=jet.colors(100), 
        axes=FALSE,
        main="")
}
dev.off()
}


### comparison of diseases for each region and sex

plot = T
if(plot){
  pdf("./results/overlap/rrho_dx_Comparisons.pdf", width = 3, height = 12)
  par(mfrow=c(4,1))
  # par(mar = c(0.2, 0.2, 0, 0))
  
  for (i in c(1:4)){
    R.obj <- RRHO(list1 = rrho.list[[i]], 
                  list2 = rrho.list[[i+4]],
                  BY = T,
                  alternative = "enrichment",
                  plots = T,
                  # log10.ind = T,
                  outputdir = "./results/overlap/rrho_dx//", 
                  labels = c(names(rrho.list)[i], names(rrho.list)[i+4])
    )
    
    image(R.obj$hypermat,
          xlab = names(rrho.list)[i],
          ylab=names(rrho.list)[i+4],
          col=jet.colors(100), 
          axes=FALSE,
          main="")
  }
  dev.off()
}

