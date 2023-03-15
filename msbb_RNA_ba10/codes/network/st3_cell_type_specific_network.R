## cell_type_specific_network.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(WGCNA,ggpubr,reshape,ComplexHeatmap,RColorBrewer,ggthemes,
                 circlize,GeneOverlap)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  load("./data/network/parameters/Combined_for_network.Rdata")
  
  ##lits of cell-type-specific RNA-markers from marker et al
  marker = read.csv("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/cell_markers.csv",skip = 2) 
  
  for(i in 1:ncol(marker)){
    marker[,i] = toupper(marker[,i])
  }
  
  #load final network
  load("./data/network/parameters/finalNetwork.RData")
  
  eigmat = MEs$eigengenes
  colnames(eigmat) = gsub("ME","M",colnames(eigmat))
  
  # Calculate Cell-Type specificity of modules
  res.p = matrix(NA, ncol = ncol(eigmat),nrow = ncol(marker))
  colnames(res.p) = colnames(eigmat)
  rownames(res.p) = colnames(marker)
  
  
  # Calculate Cell-Type specificity of modules
  res.p = matrix(NA, ncol = ncol(eigmat),nrow = ncol(marker))
  colnames(res.p) = colnames(eigmat)
  rownames(res.p) = colnames(marker)
  
  ## odd ratio results
  res.or = res.p
  
  for(cell in rownames(res.p)){
    for(mod in colnames(res.p)) {
      id1 = na.omit(gene[mods==mod])
      id2 = na.omit(unique(marker[,cell]))
      
      ## overlap
      o = newGeneOverlap(id1,id2)
      test = testGeneOverlap(o)
      # p value
      res.p[cell,mod] = signif(as.numeric(getPval(test)),2)
      # odd ratio
      res.or[cell,mod] = signif(as.numeric(getOddsRatio(test)),2)
      
    }
  }
  
  # fdr correction
  fdr = p.adjust(res.p,"fdr")
  dim(fdr) = dim(res.p); dimnames(fdr) = dimnames(res.p);
  fdr[fdr==0]=1e-250
  
  #annotation of modules
  annot_mod = data.frame(module = colnames(fdr))
  rownames(annot_mod) = colnames(fdr)
  
  # color for modules
  mod_col = names(mods)
  names(mod_col) = mods
  annot_col = list(module = mod_col)
  top_ann = HeatmapAnnotation(df = annot_mod,col = annot_col,show_legend = F,show_annotation_name = F)
  
  ## save results
  write.csv(fdr,
            file= "./results/tables/network/CellTypeEnrichmentFDR.csv")
  
  ## heatmap of cell type enrichment
  h_col = colorRamp2(breaks = c(0, max(-log10(fdr))), 
                     colors = c("white","#52006A"))
  ht_opt( legend_border = "black",
          heatmap_border = TRUE,
          annotation_border = TRUE)
  
  h1 = ComplexHeatmap::Heatmap(-log10(fdr),
                               col = h_col,
                               top_annotation = top_ann,
                               column_names_side = "top",
                               row_names_side = "left",
                               rect_gp = gpar(col = "lightgrey"),
                               # width = unit(14,"cm"),
                               # height = unit(3,"cm"),
                               show_column_names = T,
                               row_names_gp = gpar(fontsize = 8),
                               column_names_gp = gpar(fontsize = 8),
                               show_column_dend = F,
                               show_row_dend = F,
                               name = "-log10(FDR)",
                               cell_fun = function(j, i, x, y, width, height, fill) {
                                 if(fdr[i, j] < 0.05)
                                   grid.text("*", x, y,
                                             gp = gpar(fontsize = 9,
                                                       col = ifelse(-log10(fdr[i,j])>10,"white","black")))
                               })
  h1
  
  ## load correlation with covariates for plotting
  load("./data/network/mod_cov_cor.Rdata")
  corfdr = p.adjust(modulesPval)
  dim(corfdr) = dim(modulesPval);dimnames(corfdr)=dimnames(modulesPval)
  ## heatmap of correlation
  cor_col = colorRamp2(breaks = c(-0.4,0,0.4),colors = c("#DA0037","white","#4AA96C")) 
  h2 = Heatmap(moduleCors,
               name = "Correlation",
               col =cor_col,
               show_column_dend = F,
               show_row_dend = F,
               rect_gp = gpar(col = "lightgrey"),
               show_column_names = F,
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 8),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(modulesPval[i, j] < 0.05 & abs(moduleCors[i,j])>0.2)
                   grid.text("*", x, y,
                             gp = gpar(fontsize = 10,
                                       col = "black"))
               }
               )
  
  h2
  
  pdf("./results/figures/network/cell_mod_cov.pdf", 
      width = 10,height = 4)
  hlist = h1 %v% h2 
  draw(hlist,column_title = "Consensus network")
  dev.off()
  
  # automate filtering cell-specific modules
  celldat = melt(fdr)
  celldat = subset(celldat, value< 0.05)
  colnames(celldat) = c("cell","module","fdr")
  
  ## add color
  celldat$color = mod_col[match(celldat$module,names(mod_col))]
  cell_modules = celldat[order(celldat$module,decreasing = F),]
  cols = unique(celldat$color)
  names(cols) = celldat$module[match(cols,celldat$color)]
  
  saveRDS(celldat,"./codes/network/cell_modules.rds")
}

