## calculare module preservation between RNAseq modules 
## and replication datasets

if(T){
  rm(list = ls())
  library(WGCNA);library(ggplot2);library(ggthemes);library(ggrepel)
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  sex = c("female","male")
  dataset = c("TCX","CSF1","CSF2")
  res = data.frame()
  
  for (s in sex){
    for (dst in dataset){
      
      ## load discovery Mayo data
      load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
      
      discMod = mods
      discEx = exp
      discCol = colors
      rownames(discEx) = attr$gene_name 
      
      ### load replication dataset
      dir = list.files(path = "./",pattern = dst)
      load(paste0("./",dir,"/codes/network/parameters/finalNetwork_",s,".RData"))
      load(paste0("./",dir,"./codes/network/parameters/Combined_for_network.Rdata"))
      repEx = exp
      repCol = colors
      id = which(gene != "")
      gene = gene[id]
      repEx = repEx[id,]
      repCol = repCol[id]
      
      id = which(!duplicated(gene))
      gene= gene[id]
      repEx = repEx[id,]
      rownames(repEx) = gene
      repCol = repCol[id]
      
      
      setLabels = c("discovery", "rep");
      multiExpr = list(discovery = list(data = t(discEx)), 
                       rep = list(data = t(repEx)));
      multiColor = list(discovery = discCol);
      
      mp = modulePreservation(multiExpr, multiColor,
                              referenceNetworks = 1,
                              greyName = "grey",
                              dataIsExpr = T,
                              corFnc = "bicor",
                              networkType = "signed",
                              calculateQvalue = T,
                              nPermutations = 50,
                              randomSeed = 1,
                              quickCor = 0,
                              verbose = 3)
      
      ## We now analyze the data. Isolate the observed statistics and their Z scores:
      ref = 1
      test = 2
      statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
      statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
      
      #We look at the main output: the preservation medianRank and Zsummary statistics.
      # Compare preservation to quality:
      print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                  signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
      
      
      # Module labels and module sizes are also contained in the results
      modColors = rownames(mp$preservation$observed[[ref]][[test]])
      moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
      
      # Auxiliary convenience variable
      plotData = cbind(modColors,moduleSizes, mp$preservation$Z[[ref]][[test]][, 2])
      plotData = as.data.frame(plotData)
      colnames(plotData) = c("module","size","z")
      plotData = plotData[!(plotData$module %in% c("grey","gold")),]
      
      plotData$labels = discMod[match(plotData$module,names(discMod))]
      plotData$size = as.numeric(plotData$size)
      plotData$z= as.numeric(plotData$z)
      plotData = subset(plotData,z!="NaN")
      plotData = na.omit(plotData)
      
      ### save 
      plotData$sex = s
      plotData$dataset = dst
      res = rbind(res,plotData)
    }
  }
  
  ## plot results
  cols = res$module
  names(cols) = res$labels
  res$dataset = factor(res$dataset,levels = dataset)
  p = ggplot(res,aes(x = z, y = size, 
                     fill = labels),
             color = "black")+
    geom_point(shape = 21,show.legend = F,size = 3)+
    facet_grid(sex ~ dataset)+
    scale_fill_manual(values = cols)+
    scale_y_log10()+
    geom_vline(xintercept = 2,lty= 2,color = "red")+
    geom_text_repel(aes(label = labels),size =2)+
    # scale_x_log10()+
    theme_bw()+
    theme(strip.background = element_rect(fill = "darkblue"),
          strip.text = element_text(color= "white",face = "bold",size = 12))+
    labs(x = "Z summary ",
         y = "module size",
         title = "Module preservation")
  p
  ggsave(filename = "./results/figures/network/network_preservation.pdf",
         width = 8.5,height = 5,plot = p)
  
  write.csv(res, file="./results/tables/network/preservation/network_preservation.csv")
}
