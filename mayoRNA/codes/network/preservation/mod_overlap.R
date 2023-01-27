##mod_overlap.R 
if(T){
  rm(list = ls())
  
 library(GeneOverlap);library(circlize);library(ComplexHeatmap)
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  source("C:/Users/crgcomu/Desktop/Iman/Brain_meta/tools/wgnca_func.R")
  
  sex = c("female","male")
  dataset = c("TCX","CSF1","CSF2")
  # enrich = data.frame()
  match.res = data.frame()
  
  for (s in sex){
    for (dst in dataset){
      
      ## load discovery Mayo data
      load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
      mod.disc = mods
      exp.disc = exp
      attr.disc = attr
      attr.disc$mod = mod.disc
      
      ## load replication dataset
      dir = list.files(path = "./",pattern = dst)
      load(paste0("./",dir,"/codes/network/parameters/finalNetwork_",s,".RData"))
      load(paste0("./",dir,"./codes/network/parameters/Combined_for_network.Rdata"))
      
      datexp.prot = exp
      mod.prot = mods
      n.mayo = length(unique(mod.disc))
      n.prot = length(unique(mod.prot))
      
      # results
      res.p = matrix(NA, ncol = n.mayo, nrow = n.prot)
      colnames(res.p) = unique(mod.disc)
      rownames(res.p) = unique(mod.prot)
      
      # remove grey
      id = which(colnames(res.p)=="M0")
      idx = which(rownames(res.p)=="M0")
      res.p = res.p[-idx,-id]
      
      ## odd ratio results
      res.or = res.p
      
      for(col in colnames(res.p)) {
        for(row in rownames(res.p)){
          id1 = attr.disc$gene_name[col ==mod.disc]
          id2 = gene[row ==mod.prot]
          
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
      
      #annotation of modules for mayo
      annot_mayo = data.frame(module = colnames(res.p))
      rownames(annot_mayo) = colnames(res.p)
      
      # color for modules
      mayo_col = names(mod.disc)
      names(mayo_col) = mod.disc
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
      h_col = colorRamp2(breaks = c(0, 50), 
                         colors = c("white","#297F87"))
      ht_opt( legend_border = "black",
              heatmap_border = TRUE,
              annotation_border = TRUE)
      h1 = ComplexHeatmap::Heatmap(res.or,
                                   col = h_col,
                                   column_title = paste0("RNAseq- ",s),
                                   row_title = paste(dst,"proteomics- ",s,sep=" "),
                                   row_title_side = "left",
                                   # width = unit(20,"cm"),
                                   # height = unit(7,"cm"),
                                   column_names_side = "top",
                                   row_names_side = "left",
                                   top_annotation = top_ann,
                                   left_annotation = row_ann,
                                   clustering_method_columns = "mcquitty",
                                   clustering_method_rows = "mcquitty",
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
      
      pdf(paste0("./results/figures/network/preservation/RNAseq_",dst,"_",s,"_overlap.pdf"),
          width = 8, height = 3)
      draw(h1)
      dev.off()
      
      mch = matchModules(gn1 = gene, mod1 = mod.prot,
                         gn2 = attr.disc$gene_name, mod2 = mod.disc,
                         omit="grey")
      
      match_out = read.csv("./eraseMe.csv")
      colnames(match_out) = c("rep","disc","fdr")
      match_out$disc = gsub("X_","",match_out$disc)
      ## remove grey
      match_out = match_out[grep("M0",match_out$rep,invert = T),]
      ## keep top overlapped and unique
      # id = !duplicated(match_out$rep)
      # match_out = match_out[id,]
      id = !duplicated(match_out$disc)
      match_out = match_out[id,]
      match_out$sex = s
      match_out$dataset = dst
      match.res = rbind(match.res,match_out)
      file.remove("./eraseMe.csv")
      
      
    }
  }
  write.csv(match.res,file = "./results/tables/network/preservation/matched_modules_results.csv",row.names = F)
  }
  