## st3_cell_type_specific_network.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  library(ggplot2); library(biomaRt); 
  library(ggpubr); library(dplyr);library(pSI)
  library(magrittr);library(GeneOverlap);library(WGCNA)
  library(ComplexHeatmap);library(circlize)
  library(colorspace);library(reshape2)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/")
  sex = c("female","male")
  
  ##list of cell-type-specific RNA-markers from Zhang et al
  
  pSI.output <- readRDS("../mayoRNA/codes/network/zhang_pSI.rds")

    ## a df for cell-specific modules
  cell_modules = data.frame()
  
  for (s in sex){
    load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
    
    #---load final network
    eigmat = MEs$eigengenes
    colnames(eigmat) = gsub("ME","M",colnames(eigmat))
    id = which(colnames(eigmat)=="M0")
    eigmat = eigmat[,-id]
  
    res.p = matrix(NA, ncol = ncol(eigmat),nrow = 5)
    colnames(res.p) = colnames(eigmat)
    rownames(res.p) = colnames(pSI.output)
    
    for(mod in colnames(res.p)) {
      f = fisher.iteration(pSI.output, rownames(exp)[mods==mod],p.adjust = F)
      res.p[,mod] = f$`0.05 - nominal`
    }
    
    # fdr correction
    res.fdr = p.adjust(res.p,"fdr")
    dim(res.fdr) = dim(res.p); dimnames(res.fdr) = dimnames(res.p);
    res.fdr[res.fdr==0]=1e-250
    
    #annotation of modules
    annot_mod = data.frame(module = colnames(eigmat))
    rownames(annot_mod) = colnames(eigmat)
    
    # color for modules
    id = which(names(mods)== "grey")
    mod_col = unique(names(mods[-id]))
    names(mod_col) = unique(mods[-id])
    annot_col = list(module = mod_col)
    top_ann = HeatmapAnnotation(df = annot_mod,
                                col = annot_col,
                                show_legend = F,
                                show_annotation_name = F)
    
    ## save results
    out = res.fdr
    colnames(out) = paste(colnames(out), ".FDR",sep="")
    write.csv(out,
              file= paste0("./results/tables/network/CellTypeEnrichmentFDR_",s,".csv"))
    
    ## heatmap of cell type enrichment
    h_col = colorRamp2(breaks = c(0, max(-log10(res.fdr))), 
                       colors = c("white","#0C4271"))
    ht_opt( legend_border = "black",
            heatmap_border = TRUE,
            annotation_border = TRUE)
    
    h1 = ComplexHeatmap::Heatmap(-log10(res.fdr),
                                 col = h_col,
                                 top_annotation = top_ann,
                                 column_names_side = "top",
                                 row_names_side = "left",
                                 row_title = "Cell type",
                                 row_title_side = "left",
                                 width = unit(30,"cm"),
                                 height = unit(3,"cm"),
                                 show_column_names = T,
                                 row_names_gp = gpar(fontsize = 8),
                                 column_names_gp = gpar(fontsize = 8),
                                 show_column_dend = F,
                                 show_row_dend = F,
                                 name = "-log10(FDR)",
                                 cell_fun = function(j, i, x, y, width, height, fill) {
                                   if(res.fdr[i, j] < 0.001)
                                     grid.text(signif(res.fdr[i, j],1), x, y,
                                               gp = gpar(fontsize = 6, col = ifelse(-log10(res.fdr[i,j])>100,"white","black")))
                                 })
    h1
    
    ## load correlation with covariates for plotting
    cov_cor = readRDS(paste0("./codes/network/mod_cov_cor_",s,".rds"))
    rownames(cov_cor)[1] = "Condition"
    cov_p = readRDS(paste0("./codes/network/mod_cov_P_",s,".rds"))
    rownames(cov_p)[1] = "Condition"
    cor_symb = cov_p
    
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
                 width = unit(30,"cm"),
                 height = unit(4,"cm"),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(cov_p[i, j] < 0.01)
                     grid.text("*", x, y,gp = gpar(fontsize = 6))
                 })
    
    h2
    
    pdf(paste0("./results/figures/network/cell_mod_cov_",s,".pdf"), width = 14,height = 4)
    hlist = h1 %v% h2 
    draw(hlist,column_title =paste("RNAseq Mayo-",s, " Network",sep=""))
    dev.off()
    
    # automate filtering cell-specific modules
    celldat = melt(res.fdr)
    celldat = subset(celldat,value <1e-10)
    colnames(celldat) = c("cell","module","fdr")
    celldat$sex = s
    celldat$color = mod_col[match(celldat$module,names(mod_col))]
    cell_modules = rbind(cell_modules,celldat)
    cell_modules = cell_modules[order(cell_modules$module,decreasing = F),]
    cols = unique(cell_modules$color)
    names(cols) = cell_modules$module[match(cols,cell_modules$color)]
    
  }
  saveRDS(cell_modules,"./codes/network/cell_modules.rds")
}
