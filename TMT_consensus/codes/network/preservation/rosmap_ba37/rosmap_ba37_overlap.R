##mod_overlap.R 
if(T){
  rm(list = ls())
  
  pacman::p_load(GeneOverlap,circlize,ComplexHeatmap)
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
  source("C:/Users/crgcomu/Desktop/Iman/Brain_meta/tools/wgnca_func.R")
  
  match.res = data.frame()
      
      ## load discovery TMT rnaein
      load("./data/network/parameters/finalNetwork.RData")
      mod.disc = mods
      exp.disc = datExp
      attr.disc = data.frame(gene = rownames(datExp),mod.disc)
      
      ## load replication dataset
      load("../Rosmap_ba37_ba6/data/network/parameters/finalNetwork.RData")
      load("../Rosmap_ba37_ba6/data/normalized.Rdata")
      
      datexp.tcx = datExp
      mod.tcx = mods
      attr.tcx = data.frame(gene = rownames(datExp),
                            mod = mod.tcx)
      n.prot = length(unique(mod.disc))
      n.tcx = length(unique(mod.tcx))
      
      # results
      res.p = matrix(NA, ncol = n.prot, nrow = n.tcx)
      colnames(res.p) = unique(mod.disc)
      rownames(res.p) = unique(mod.tcx)
      
      # remove grey
      id = which(colnames(res.p)!="M0")
      idx = which(rownames(res.p)!="M0")
      res.p = res.p[idx,id]
      
      ## odd ratio results
      res.or = res.p
      
      for(col in colnames(res.p)) {
        for(row in rownames(res.p)){
          id1 = attr.disc$gene[col == mod.disc]
          id2 = attr.tcx$gene[row == mod.tcx]
          
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
      res.fdr = p.adjust(res.p,"fdr")
      dim(res.fdr) = dim(res.p); dimnames(res.fdr) = dimnames(res.p);
      res.fdr[res.fdr ==0] = 1e-250
      sym = res.fdr
      sym[sym<0.05]="*"
      sym[sym<0.01]="**"
      sym[sym<0.001 ]="***"
      sym[sym>0.05]=""
      
      #annotation of modules for prot
      annot_prot = data.frame(module = colnames(res.p))
      rownames(annot_prot) = colnames(res.p)
      
      # color for modules
      prot_col = names(mod.disc)
      names(prot_col) = mod.disc
      annot_col_prot = list(module = prot_col)
      top_ann = HeatmapAnnotation(df = annot_prot,col = annot_col_prot,
                                  show_legend = F,show_annotation_name = F)
      
      
      #annotation of modules for rna
      annot_rna = data.frame(module = rownames(res.p))
      rownames(annot_rna) = rownames(res.p)
      
      # color for modules
      rna_col = names(mod.tcx)
      names(rna_col) = mod.tcx
      annot_col_rna = list(module = rna_col)
      row_ann = rowAnnotation(df = annot_rna,col = annot_col_rna,
                              show_legend = F,show_annotation_name = F)
      
      ## heatmap of cell type enrichment
      h_col = colorRamp2(breaks = c(0, 50), 
                         colors = c("white","#297F87"))
      ht_opt( legend_border = "black",
              heatmap_border = TRUE,
              annotation_border = TRUE)
      h1 = ComplexHeatmap::Heatmap(res.or,
                                   col = h_col,
                                   column_title = "TMT consensus proteomics",
                                   row_title = "Temporal cortex proteomics",
                                   row_title_side = "left",
                                   # width = unit(20,"cm"),
                                   # height = unit(7,"cm"),
                                   column_names_side = "top",
                                   row_names_side = "left",
                                   top_annotation = top_ann,
                                   left_annotation = row_ann,
                                   # clustering_method_columns = "mcquitty",
                                   # clustering_method_rows = "mcquitty",
                                   row_names_gp = gpar(fontsize = 7),
                                   column_names_gp = gpar(fontsize = 7),
                                   cluster_columns = T,
                                   show_column_dend = F,
                                   show_row_dend = F,
                                   rect_gp = gpar(col = "#F3F1F5"),
                                   name = "OR",
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                     if(res.fdr[i, j] < 0.05)
                                       grid.text(sym[i, j], x, y,gp = gpar(fontsize = 8))
                                   })
      
      pdf("./results/figures/network/preservation/rosmap_ba37_overlap.pdf",
          width = 8, height = 6)
      draw(h1)
      dev.off()
      
      mch = matchModules(gn1 = attr.tcx$gene, mod1 = mod.tcx,
                         gn2 = attr.disc$gene, mod2 = mod.disc,
                         omit="grey")
      
      match_out = read.csv("./eraseMe.csv")
      colnames(match_out) = c("Tcx","TMT","fdr")
      match_out$TMT = gsub("X_","",match_out$TMT)
      
      ## remove grey
      match_out = match_out[grep("M0",match_out$Tcx,invert = T),]
      
      ## keep top overlapped and unique
      # id = !duplicated(match_out$rep)
      # match_out = match_out[id,]
      # id = !duplicated(match_out$Protein)
      # match_out = match_out[id,]
      file.remove("./eraseMe.csv")
      ggplot(match_out,aes(x = Tcx,y = TMT,fill = -log10(fdr)))+
        geom_tile()+
        scale_fill_continuous(type = "viridis")+
        theme_bw()+
        theme(axis.text.x.bottom = element_text(angle = 45))
      
  write.csv(match_out,file = "./results/tables/network/preservation/rosmap_BA37_matched_modules.csv",row.names = F)
}
