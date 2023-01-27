## cell_type_specific_network.R
if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  library(ggplot2); library(biomaRt); 
  library(ggpubr); library(dplyr);library(reshape2)
  library(magrittr);library(GeneOverlap);
  library(ComplexHeatmap);library(circlize)
  library(colorspace)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/replication_proteomics_CSF2")
  sex = c("female","male")
  ## a df for cell-specific modules
  cell_modules = data.frame()
  
  ##lits of cell-type-specific RNA-markers from Zhang et al
  zhang = read.csv("../../data/scRNAseq/zhang/list.csv",skip=2) 
  zhang = zhang[,8:11]
  colnames(zhang) = c("astrocyte","micro","neuron","oligo")
  
  for(i in 1:ncol(zhang)){
    zhang[,i] = toupper(zhang[,i])
  }
  
  for (s in sex){
    load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
    
    #---load final network
    
    eigmat = MEs$eigengenes
    colnames(eigmat) = gsub("ME","M",colnames(eigmat))
    id = which(colnames(eigmat)=="M0")
    eigmat = eigmat[,-id]
    
    # Calculate Cell-Type specificity of modules
    res.p = matrix(NA, ncol = ncol(eigmat),nrow = ncol(zhang))
    colnames(res.p) = colnames(eigmat)
    rownames(res.p) = colnames(zhang)
    
    ## odd ratio results
    res.or = res.p
    
    for(mod in colnames(res.p)) {
      for(cell in rownames(res.p)){
        id1 = gene[mods==mod]
        print(id1)
        id2 = unique(zhang[,cell])
        
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
    out = fdr
    colnames(out) = paste(colnames(out), ".FDR",sep="")
    write.csv(out,
              file= paste0("./results/tables/network/CellTypeEnrichmentFDR_",s,".csv"))
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
                                 width = unit(15,"cm"),
                                 height = unit(3,"cm"),
                                 show_column_names = T,
                                 row_names_gp = gpar(fontsize = 8),
                                 column_names_gp = gpar(fontsize = 8),
                                 show_column_dend = F,
                                 show_row_dend = F,
                                 name = "-log10(FDR)",
                                 cell_fun = function(j, i, x, y, width, height, fill) {
                                   if(fdr[i, j] < 0.05)
                                     grid.text(signif(fdr[i, j],1), x, y,
                                               gp = gpar(fontsize = 8,
                                                         col = ifelse(-log10(fdr[i,j])>20,"white","black")))
                                 })
    h1
    
    ## load correlation with covariates for plotting
    cov_cor = readRDS(paste0("./codes/network/mod_cov_cor_",s,".rds"))
    rownames(cov_cor)[1] = "Condition"
    cov_p = readRDS(paste0("./codes/network/mod_cov_P_",s,".rds"))
    rownames(cov_p)[1] = "Condition"
    cor_symb = cov_p
    
    ## heatmap of correlation
    cor_col = colorRamp2(breaks = c(-0.3,0,0.3),colors = c("#DA0037","white","#4AA96C")) 
    h2 = Heatmap(cov_cor,
                 name = "Correlation",
                 col =cor_col,
                 show_column_dend = F,
                 show_row_dend = F,
                 rect_gp = gpar(col = "lightgrey"),
                 show_column_names = F,
                 row_names_side = "left",
                 row_names_gp = gpar(fontsize = 8),
                 width = unit(15,"cm"),
                 height = unit(4,"cm"),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(cov_p[i, j] < 0.01)
                     grid.text(signif(cov_p[i, j],1), x, y,
                               gp = gpar(fontsize = 8,
                                         col = ifelse(-log10(cov_p[i,j])>30,"white","black")))
                 })
    
    h2
    
    pdf(paste0("./results/figures/network/cell_mod_cov_",s,".pdf"), 
        width = 8,height = 4)
    hlist = h1 %v% h2 
    draw(hlist,column_title =paste("CSF2 -",s, " network",sep=""))
    dev.off()
    
    # automate filtering cell-specific modules
    celldat = melt(fdr)
    celldat = subset(celldat, value< 0.05)
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

# save.image("./codes/network/cell_type_enrich.Rdata")
