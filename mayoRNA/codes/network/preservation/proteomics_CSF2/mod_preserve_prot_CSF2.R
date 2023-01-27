##mod_preserve_mayo_proteomics.R 
## calculare module preservation between discovery and replication RNAseq Temporal-MSBB
if(T){
  rm(list = ls())
  library(WGCNA)
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  sex = c("female","male")
  res = data.frame()
  
  for (s in sex){
    ## load discovery Mayo data
    load(paste0("./codes/network/finalNetwork_",s,".RData",sep=""))
    discMod = mods
    discEx = exp
    discCol = colors
    rownames(discEx) = attr$gene_name 
    
    ## load Mayo proteomis network 
    load(paste0("./replication_proteomics_CSF2/codes/network/parameters/finalNetwork_",s,".RData"))
    repEx = exp
    repCol = colors
    # id = which(gene == "")
    # gene = gene[-id]
    # repEx = repEx[-id,]
    # repCol = repCol[-id]
    
    # id = which(duplicated(gene))
    # gene= gene[-id]
    # repEx = repEx[-id,]
    # rownames(repEx) = gene
    # repCol = repCol[-id]
    
    
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
    
    # Save the results
    # save(mp, file = "./codes/network/preservation/mayoProteomics_modulePreservation.RData");
    
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
    plotData = subset(plotData,z!="Na")
    plotData = na.omit(plotData)
    
    ## ploting
    pdf(paste0("./results/figures/network/preservation/CSF2_preservation_",s,".pdf"),
        width = 6, height = 4)
    plot(x = plotData[,"z"],
         y = plotData[,"size"],
         col = 1, 
         bg = plotData$module,
         pch = 21,
         main = paste0("CSF2 proteomics- ",s),
         cex = 2.4,
         xlab =  expression("Z "[summary]),
         ylab = "Module size",
         log = "y",
         xlim = c(0,max(plotData[,"z"])),
         ylim = c(2, max(plotData[,"size"])), 
         cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
    labelPoints(x = plotData$z, y = plotData$size, labels = plotData$labels, cex = 0.7, offs = 0.08);
    # For Zsummary, add threshold lines
    abline(v=0)
    abline(v=2, col = "blue", lty = 5)
    abline(v=10, col = "red", lty = 5)
    dev.off()
    
    ### save 
    plotData$sex = s
    plotData$dataset = "CSF2"
    res = rbind(res,plotData)
  }
  saveRDS(res,file = paste("./results/tables/network/preservation/csf2_preserve",s,".rds"))
  
}
