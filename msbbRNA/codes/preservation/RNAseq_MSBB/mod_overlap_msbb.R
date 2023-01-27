##mod_overlap_msbb.R 
library(WGCNA);library(GeneOverlap);library(circlize);library(ComplexHeatmap)

## load discovery Mayo data
load("./codes/network/parameters/finalNetwork.RData")
mods.discovery = mods
datexp.discovery = datExp

## load MSBB data
load("./replication_msbb/codes/network/parameters/finalNetwork.RData")
datexp.msbb = datExp
mod.msbb = mods
n.mayo = length(unique(mods.discovery))
n.msbb = length(unique(mod.msbb))

# results
res.p = matrix(NA, ncol = n.mayo, nrow = n.msbb)
colnames(res.p) = unique(mods.discovery)
rownames(res.p) = unique(mod.msbb)

# remove grey
res.p = res.p[-which(rownames(res.p)=="M0"),-which(colnames(res.p)=="M0")]

## odd ratio results
res.or = res.p

for(col in colnames(res.p)) {
  for(row in rownames(res.p)){
    id1 = rownames(datexp.discovery)[col ==mods.discovery]
    id2 = rownames(datexp.msbb)[row ==mod.msbb]
    
    ## overlap
    o = newGeneOverlap(id1,id2)
    test = testGeneOverlap(o)
    # p value
    res.p[row,col] = signif(as.numeric(getPval(test)),2)
    # odd ratio
    res.or[row,col] = signif(as.numeric(getOddsRatio(test)),2)
    
  }
}

# fdr correction
res.p.fdr = p.adjust(res.p,"fdr")
dim(res.p.fdr) = dim(res.p); dimnames(res.p.fdr) = dimnames(res.p);
res.p.fdr[res.p.fdr < 1e-300] = 1e-300
sym = res.p.fdr
sym[sym<0.05]="*"
sym[sym<0.01]="**"
sym[sym<0.001]="***"

## save results
out = res.p.fdr
colnames(out) = paste(colnames(out), ".FDR",sep="")
write.csv(out,file="./results/tables/network/CellTypeEnrichmentFDR.csv")

#annotation of modules for mayo
annot_mayo = data.frame(module = colnames(res.p))
rownames(annot_mayo) = colnames(res.p)

# color for modules
id = which(names(mods.discovery)=="grey")
mayo_col = unique(names(mods.discovery[-id]))
names(mayo_col) = unique(mods.discovery[-id])
annot_col_mayo = list(module = mayo_col)
top_ann = HeatmapAnnotation(df = annot_mayo,col = annot_col_mayo,
                            show_legend = F,show_annotation_name = F)

#annotation of modules for msbb
annot_msbb = data.frame(module = rownames(res.p))
rownames(annot_msbb) = rownames(res.p)

# color for modules
id = which(names(mod.msbb)=="grey")
msbb_col = unique(names(mod.msbb[-id]))
names(msbb_col) = unique(mod.msbb[-id])
annot_col_msbb = list(module = msbb_col)
row_ann = rowAnnotation(df = annot_msbb,col = annot_col_msbb,
                        show_legend = F,show_annotation_name = F)

## heatmap of cell type enrichment
h_col = colorRamp2(breaks = c(0, max(-log10(res.p.fdr))), 
                   colors = c("white","purple"))
ht_opt( legend_border = "black",
        heatmap_border = TRUE,
        annotation_border = TRUE)

h1 = ComplexHeatmap::Heatmap(-log10(res.p.fdr),
                             col = h_col,
                             column_title = expression("Discovery "[mayo]),
                             row_title = expression("Replication "[msbb]),
                             row_title_side = "left",
                             width = unit(20,"cm"),
                             height = unit(7,"cm"),
                             column_names_side = "top",
                             row_names_side = "left",
                             top_annotation = top_ann,
                             left_annotation = row_ann,
                             clustering_method_columns = "ward.D2",
                             clustering_method_rows = "ward.D2",
                             row_names_gp = gpar(fontsize = 6),
                             column_names_gp = gpar(fontsize = 6),
                             border = "black",
                             cluster_columns = T,
                             show_column_dend = F,
                             show_row_dend = F,
                             # rect_gp = gpar(col = "lightgrey"),
                             name = "-log10(FDR)",
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(res.p.fdr[i, j] < 0.05)
                                 grid.text(sym[i, j], x, y,gp = gpar(fontsize = 7))
                             })
pdf("./results/figures/network/preservation/mayo_msbb_network_overlap.pdf",
    width = 10, height = 4)
draw(h1)
dev.off()
