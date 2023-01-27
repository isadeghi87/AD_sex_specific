if(T){
  rm(list = ls())
  library(WGCNA);library(GeneOverlap);library(circlize);library(ComplexHeatmap)
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  sex = c("female","male")
  enrich = data.frame()
  res.list = list()
  
  for (s in sex){
    ## load discovery Mayo data
    load(paste0("./codes/network/finalNetwork_",s,".RData",sep=""))
    mods.discovery = mods
    datexp.discovery = exp
    attr.disc = attr
    
    ## load CSF2 proteomics data
    load(paste0("./replication_proteomics_CSF2/codes/network/parameters/finalNetwork_",s,".RData"))
    load("./replication_proteomics_CSF2/codes/network/parameters/Combined_for_network.Rdata")
    
    datexp.prot = exp
    mod.prot = mods
    n.mayo = length(unique(mods.discovery))
    n.prot = length(unique(mod.prot))
    
    # results
    res.p = matrix(NA, ncol = n.mayo, nrow = n.prot)
    colnames(res.p) = unique(mods.discovery)
    rownames(res.p) = unique(mod.prot)
    
    # remove grey
    id = which(colnames(res.p)=="M0")
    idx = which(rownames(res.p)=="M0")
    res.p = res.p[-idx,-id]
    
    ## odd ratio results
    res.or = res.p
    
    for(col in colnames(res.p)) {
      for(row in rownames(res.p)){
        id1 = attr.disc$gene_name[col ==mods.discovery]
        id2 = gene[row ==mod.prot]
        size = length(union(attr$gene_name,gene))
        ## overlap
        o = newGeneOverlap(id1,id2,
                           # spec = "hg19.gene",
                           genome.size = size
                           )
        test = testGeneOverlap(o)
        # p value
        res.p[row,col] = signif(as.numeric(getPval(test)),2)
        # odd ratio
        res.or[row,col] = signif(as.numeric(getOddsRatio(test)),2)
        
      }
    }
    
    # fdr correction
    res.fdr = p.adjust(res.p,"fdr")
    dim(res.fdr) = dim(res.p); dimnames(res.fdr) = dimnames(res.p);
    res.fdr[res.fdr ==0] = 1e-250
    res.p[res.p ==0] = 1e-250
    sym = res.p
    sym[sym<0.05 & res.or>2]="*"
    sym[sym<0.01 & res.or>5]="**"
    sym[sym<0.001 & res.or>10]="***"
    sym[sym>0.05 | res.or<=2]=""
    
    ## save results
    res.list[["fdr"]] = res.fdr
    res.list[["or"]] = res.p
    res.list[["sym"]] = res.or
    
    saveRDS(res.list,file= paste0("./results/tables/network/preservation/mayoRNA_csf2_overlap_",s,".rds"))
    
    #annotation of modules for mayo
    annot_mayo = data.frame(module = colnames(res.p))
    rownames(annot_mayo) = colnames(res.p)
    
    # color for modules
    mayo_col = names(mods.discovery)
    names(mayo_col) = mods.discovery
    annot_col_mayo = list(module = mayo_col)
    top_ann = HeatmapAnnotation(df = annot_mayo,col = annot_col_mayo,
                                show_legend = F,show_annotation_name = F)
    
    
    #annotation of modules for prot
    annot_prot = data.frame(module = rownames(res.p))
    rownames(annot_prot) = rownames(res.p)
    
    # color for modules
    prot_col = names(mod.prot)
    names(prot_col) = mod.prot
    annot_col_prot = list(module = prot_col)
    row_ann = rowAnnotation(df = annot_prot,col = annot_col_prot,
                            show_legend = F,show_annotation_name = F)
    
    ## heatmap of cell type enrichment
    h_col = colorRamp2(breaks = c(0, max(-log10(res.p))), 
                       colors = c("white","#297F87"))
    ht_opt( legend_border = "black",
            heatmap_border = TRUE,
            annotation_border = TRUE)
    coltitle = expr()
    h1 = ComplexHeatmap::Heatmap(-log10(res.p),
                                 col = h_col,
                                 column_title = paste0("Mayo RNAseq- ",s),
                                 row_title = paste0("CSF2 proteomics- ",s),
                                 row_title_side = "left",
                                 width = unit(16,"cm"),
                                 height = unit(7,"cm"),
                                 column_names_side = "top",
                                 row_names_side = "left",
                                 top_annotation = top_ann,
                                 left_annotation = row_ann,
                                 # clustering_method_columns = "ward.D2",
                                 # clustering_method_rows = "ward.D2",
                                 row_names_gp = gpar(fontsize = 7),
                                 column_names_gp = gpar(fontsize = 7),
                                 cluster_columns = T,
                                 show_column_dend = F,
                                 show_row_dend = F,
                                 rect_gp = gpar(col = "#F3F1F5"),
                                 name = "-log10(p)",
                                 cell_fun = function(j, i, x, y, width, height, fill) {
                                   if(res.p[i, j] < 0.05)
                                     grid.text(sym[i, j], x, y,gp = gpar(fontsize = 7))
                                 })
    
    pdf(paste0("./results/figures/network/preservation/mayoRNA_CSF2_overlap_",s,".pdf"),
        width = 8, height = 4)
    draw(h1)
    dev.off()
    
    df = melt(res.or)
    p = melt(res.fdr)
    df$p = p$value
    df$Var1 = as.character(df$Var1)
    df = df[df$value>2 & df$p < 0.05,]
    df$sex = s
    df$dataset = "CSF2-proteomics"
    colnames(df)[1:3] = c("discovery","rep","OR")
    enrich = rbind(enrich,df) 
    
  }
  saveRDS(enrich,"./results/tables/network/preservation/csf2_overlap.rds")
}