## cell_type_specific_network.R
rm(list=ls()); options(stringsAsFactors = F)
library(ggplot2); library(biomaRt); 
library(ggpubr); library(dplyr);
library(magrittr);library(GeneOverlap);
library(ComplexHeatmap);library(circlize)
library(colorspace)
condition = TRUE

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/bulk_analysis/replication_msbb/")

#---load final network
load("./codes/network/parameters/finalNetwork.RData")

eigmat = MEs$eigengenes
colnames(eigmat) = gsub("ME","M",colnames(eigmat))
eigmat = eigmat[,-1]
kME = signedKME(t(datExp), MEs$eigengenes)
colnames(kME) = gsub("kME", "M", colnames(kME))
kME = kME[,-1]

##lits of cell-type-specific RNA-markers from Zhang et al
zhang = read.csv("../../data/scRNAseq/zhang/list.csv",skip=2) 
zhang = zhang[,2:6]
colnames(zhang)[3] = "Oligo"

for(i in 1:ncol(zhang)){
  zhang[,i] = toupper(zhang[,i])
}

#----Annotating gene IDs
getinfo <- c("ensembl_gene_id",
             "external_gene_name")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
attr <- getBM(
  attributes = getinfo,
  filters = "ensembl_gene_id",
  values = rownames(datExp),
  mart = mart)

# Calculate Cell-Type specificity of modules
res.p = matrix(NA, ncol = ncol(eigmat),nrow = 5)
colnames(res.p) = colnames(eigmat)
rownames(res.p) = colnames(zhang)

## odd ratio results
res.or = res.p

for(mod in colnames(res.p)) {
  for(cell in rownames(res.p)){
    id1 = attr$external_gene_name[mod ==mods]
    id2 = zhang[,cell]
    id2 = toupper(id2[id2!= ""])
    
    ## overlap
    o = newGeneOverlap(id1,id2,spec = "hg19.gene")
    test = testGeneOverlap(o)
    # p value
    res.p[cell,mod] = signif(as.numeric(getPval(test)),2)
    # odd ratio
    res.or[cell,mod] = signif(as.numeric(getOddsRatio(test)),2)
    
  }
}

# fdr correction
res.p.fdr = p.adjust(res.p,"fdr")
dim(res.p.fdr) = dim(res.p); dimnames(res.p.fdr) = dimnames(res.p);

#annotation of modules
annot_mod = data.frame(module = colnames(eigmat))
rownames(annot_mod) = colnames(eigmat)

# color for modules
mod_col = unique(names(mods))[-2]
names(mod_col) = unique(mods)[-2]
annot_col = list(module = mod_col)
top_ann = HeatmapAnnotation(df = annot_mod,col = annot_col,show_legend = F,show_annotation_name = F)

## save results
out = res.p.fdr
colnames(out) = paste(colnames(out), ".FDR",sep="")
write.csv(out,file="./results/tables/network/CellTypeEnrichmentFDR.csv")

## heatmap of cell type enrichment
h_col = colorRamp2(breaks = c(0, max(-log10(res.p.fdr))), 
                   colors = c("white","blue"))
ht_opt( legend_border = "black",
        heatmap_border = TRUE,
        annotation_border = TRUE)

h1 = ComplexHeatmap::Heatmap(-log10(res.p.fdr),
                             col = h_col,
                             top_annotation = top_ann,
                             column_names_side = "top",
                             row_names_side = "left",
                             row_title = "Cell type",
                             row_title_side = "left",
                             width = unit(15,"cm"),
                             height = unit(3,"cm"),
                             show_column_names = T,
                             row_names_gp = gpar(fontsize = 8),
                             column_names_gp = gpar(fontsize = 8),
                             show_column_dend = F,
                             show_row_dend = F,
                             name = "-log10(FDR)",
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(res.p.fdr[i, j] < 0.001)
                                 grid.text(signif(res.p.fdr[i, j],1), x, y,gp = gpar(fontsize = 7))
                             })
h1

## load correlation with covariates for plotting
cov_cor = readRDS("./codes/network/mod_cov_cor.rds")
rownames(cov_cor)[1] = "Condition"
cov_p = readRDS("./codes/network/mod_cov_P.rds")
rownames(cov_p)[1] = "Condition"

## heatmap of correlation
cor_col = colorRamp2(range(cov_cor),c("white","orange")) 
h2 = Heatmap(cov_cor,
             name = "Correlation",
             col =cor_col,
             show_column_dend = F,
             show_row_dend = F,
             row_title = "Covariates",
             row_title_side = "left",
             show_column_names = F,
             row_names_side = "left",
             row_names_gp = gpar(fontsize = 8),
             width = unit(14,"cm"),
             height = unit(4,"cm"),
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(cov_p[i, j] < 0.01)
                 grid.text(signif(cov_p[i, j],1), x, y,gp = gpar(fontsize = 7))
             })

h2

pdf("./results/figures/network/cell_mod_cov.pdf", width = 8,height = 5)
hlist = h1 %v% h2 
draw(hlist,column_title = expression("MSBB "[replication]*" network"))
dev.off()

save.image("./codes/network/cell_type_enrich.Rdata")
