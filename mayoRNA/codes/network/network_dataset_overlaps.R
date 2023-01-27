if(T){
  rm(list=ls())
  library(ggplot2)
  library(RColorBrewer)
  library(ggthemes)
  library(circlize)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  sex = c("female","male")
  
  res= read.csv("./results/tables/network/preservation/matched_modules_results.csv")
  
  for(s in sex){
    load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
    
    df = res[res$sex==s,]
    df = df[,c("disc","dataset","fdr")]
    df = pivot_wider(df,names_from =dataset,values_from  =fdr)
    mat = t(df[,2:4])
    colnames(mat) = df$disc
    mat[is.na(mat)]=1
    mat[mat==0]=1e-300
    
    annot_mod = data.frame(RNAseq = colnames(mat))
    rownames(annot_mod) = colnames(mat)
    # color for modules
    mod_col = unique(names(mods))
    names(mod_col) = unique(mods)
    mod_col = mod_col[match(colnames(mat),names(mod_col))]
    annot_col = list(RNAseq = mod_col)
    top_ann = HeatmapAnnotation(df = annot_mod, col = annot_col,show_legend = F,show_annotation_name = F)
    
    
    col_fun = colorRamp2(c(0, -log10(1e-300)), 
                         c( "white", ifelse(s=="male","blue","red")))
    h = Heatmap(-log10(mat),
            col = col_fun,
            top_annotation = top_ann,
            column_names_side = "top",
            row_names_side = "left",
            cluster_columns = F,cluster_rows = F,
            column_title = s,
            rect_gp = gpar(col = "lightgrey"),
            show_column_names = T,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8),
            show_column_dend = F,
            show_row_dend = F,
            name = "-log10(FDR)",
            cell_fun = function(j, i, x, y, width, height, fill) {
              if(mat[i, j] < 0.05)
                grid.text("*", x, y,
                          gp = gpar(fontsize =10,
                                    col = ifelse(-log10(mat[i,j])>100,"white","black")))
            })
    pdf(file = paste0("./results/figures/network/network_datasets_overlap_",s,".pdf"),
        width = 5,height = 2)
    draw(h)
    dev.off()
    
  }
 
}