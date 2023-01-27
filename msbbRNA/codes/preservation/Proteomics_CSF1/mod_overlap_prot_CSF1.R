##mod_overlap_mprot_CSF1.R
## overlap of module networks from Mayo RNAseq and CSF1 proteomics

library(WGCNA);library(GeneOverlap);library(circlize);library(ComplexHeatmap)
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/bulk_analysis")

## load discovery Mayo RNA data
load("./codes/network/parameters/finalNetwork.RData")
mods.discovery = mods
datexp.discovery = datExp

## load CSF1 proteomics data
load("./replication_proteomics_CSF1/codes/network/parameters/finalNetwork.RData")
datexp.prot = datExp
mod.prot = mods
n.mayo = length(unique(mods.discovery))
n.prot = length(unique(mod.prot))

# results
res.p = matrix(NA, ncol = n.mayo, nrow = n.prot)
colnames(res.p) = unique(mods.discovery)
rownames(res.p) = unique(mod.prot)

# remove grey
res.p = res.p[which(rownames(res.p) != "M0"), which(colnames(res.p) != "M0")]

## odd ratio results
res.or = res.p

for(col in colnames(res.p)) {
  for(row in rownames(res.p)){
    id1 = rownames(datexp.discovery)[col ==mods.discovery]
    id2 = rownames(datexp.prot)[row ==mod.prot]
    
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
# res.p.fdr = res.p
# dim(res.p.fdr) = dim(res.p); dimnames(res.p.fdr) = dimnames(res.p);
# res.p.fdr[res.p.fdr < 1e-300] = 1e-300
sym = res.p
sym[sym<0.05]="*"
sym[sym<0.01]="**"
sym[sym<0.001]="***"

## save results
write.csv(res.p,file="./results/tables/network/preservation/mayoRNA_csf1Protein_modOverlap.csv")

#annotation of modules for mayo
annot_mayo = data.frame(module = colnames(res.p))
rownames(annot_mayo) = colnames(res.p)

# color for modules
id = which(unique(mods.discovery) != "M0")
mayo_col = unique(names(mods.discovery))[id]
names(mayo_col) = unique(mods.discovery)[id]
annot_col_mayo = list(module = mayo_col)
top_ann = HeatmapAnnotation(df = annot_mayo,col = annot_col_mayo,
                            show_legend = F,show_annotation_name = F)

#annotation of modules for prot
annot_prot = data.frame(module = rownames(res.p))
rownames(annot_prot) = rownames(res.p)

# color for modules
id = which(unique(mod.prot) != "M0")
prot_col = unique(names(mod.prot))[id]
names(prot_col) = unique(mod.prot)[id]
annot_col_prot = list(module = prot_col)
row_ann = rowAnnotation(df = annot_prot,col = annot_col_prot,
                        show_legend = F,show_annotation_name = F)

## heatmap of cell type enrichment
h_col = colorRamp2(breaks = c(0, max(-log10(res.p))), 
                   colors = c("white","#3E00FF"))
ht_opt( legend_border = "black",
        heatmap_border = TRUE,
        annotation_border = TRUE)

h1 = ComplexHeatmap::Heatmap(-log10(res.p),
                             col = h_col,
                             column_title = expression("RNAseq"[mayo]*" Network"),
                             row_title = expression("CSF Proteomics Network"),
                             row_title_side = "left",
                             width = unit(20,"cm"),
                             height = unit(7,"cm"),
                             column_names_side = "top",
                             row_names_side = "left",
                             top_annotation = top_ann,
                             left_annotation = row_ann,
                             clustering_method_columns = "ward.D2",
                             clustering_method_rows = "ward.D2",
                             row_names_gp = gpar(fontsize = 7),
                             column_names_gp = gpar(fontsize = 7),
                             cluster_columns = T,
                             show_column_dend = F,
                             show_row_dend = F,
                             rect_gp = gpar(col = "#F3F1F5"),
                             name = "-log10(FDR)",
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(res.p[i, j] < 0.05)
                                 grid.text(sym[i, j], x, y,
                                           gp = gpar(fontsize = 8,
                                                     col = ifelse(-log10(res.p[i,j])>4,"white","black")))
                             })

h1
pdf("./results/figures/network/preservation/mayoRNA_csfProtein_modOverlap.pdf",
    width = 10, height = 4)
draw(h1)
dev.off()
