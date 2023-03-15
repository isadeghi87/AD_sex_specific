## calculare module preservation between ROSMAP protein and RNA
if(T){
  rm(list = ls())
  pacman::p_load(WGCNA,ggplot2)
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
    
  ## load TMT proteomics
    load("./data/network/parameters/finalNetwork.RData")
    discMod = mods
    discEx = datExp
    discCol = colors
    
    ## load RNA network 
    load("../Rosmap_ba37_ba6/data/network/parameters/finalNetwork.RData")
    repEx = datExp
    repCol = colors
    
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
    plotData$preserved = ifelse(plotData$z>1.96,"TRUE","FALSE")
    table(plotData$preserved)
    
    ## ploting
    cols = plotData$module
    names(cols) = plotData$labels
    library(ggrepel)
    p = ggplot(plotData,aes(x = size, y = z, 
                       fill = labels),
               color = "black")+
      geom_point(shape = 21,show.legend = F,size = 4)+
      scale_fill_manual(values = cols)+
      scale_x_log10()+
      geom_hline(yintercept = 1.96,lty= 2,color = "red")+
      geom_text_repel(aes(label = labels),size =3)+
      labs(y = "Z summary ",
           x = "module size",
           title = "Module preservation in BA37 proteomics")+
      theme_bw()+
      theme(strip.background = element_rect(fill = "darkblue"),
            plot.title = element_text(hjust = 0.5),
            strip.text = element_text(color= "white",face = "bold",size = 12))
     
    p
    ggsave("./results/figures/network/rosmap_ba37_preserved.pdf",
           plot = p,
           width = 6,height =  4)
    
  ### save 
  write.csv(x = plotData,"./results/tables/network/preservation/rosmap_BA37.csv")
}
