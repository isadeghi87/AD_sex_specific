##mod_preserve_msbb.R 
## calculare module preservation between discovery and replication RNAseq Temporal-MSBB
library(WGCNA)
load("./codes/network/parameters/finalNetwork.RData")
mods.discovery = mods
datexp.discovery = datExp
colors.discovery = colors

## load MSBB data
load("./replication_msbb/codes/network/parameters/finalNetwork.RData")
datexp.msbb = datExp
colors.msbb = colors

setLabels = c("discovery", "msbb");
multiExpr = list(discovery = list(data = t(datexp.discovery)), 
                 msbb = list(data = t(datexp.msbb)));
multiColor = list(discovery = colors.discovery);

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,greyName = "grey",
                          dataIsExpr = T,corFnc = "bicor",
                          networkType = "signed",
                          calculateQvalue = T,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );

# Save the results
save(mp, file = "./codes/network/preservation/modulePreservation.RData");

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

plotData = cbind(modColors,moduleSizes, mp$preservation$Z[[ref]][[test]][, 2])
plotData = as.data.frame(plotData)
colnames(plotData) = c("module","size","z")
plotData = plotData[!(plotData$module %in% c("grey","gold")),]

plotData$labels = mods.discovery[match(plotData$module,names(mods.discovery))]
plotData$size = as.numeric(plotData$size)
plotData$z= as.numeric(plotData$z)

# Auxiliary convenience variable

## ploting
if(T){
  pdf("./results/figures/network/preservation/Mayo_MSBB_preservation.pdf",
      width = 8, height = 6)
plot(x = plotData[,"z"],
     y = plotData[,"size"],
     col = 1, 
     bg = plotData$module,
     pch = 21,
       main = expression("Module Preservataion"),
       cex = 2.4,
       xlab =  expression("Z "[summary]),
     ylab = "Module size",
     log = "y",
       xlim = c(0,max(plotData[,"z"])),
       ylim = c(50, 1000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(x = plotData$z, y = plotData$size, labels = plotData$labels, cex = 0.8, offs = 0.08);
  # For Zsummary, add threshold lines
    abline(v=0)
    abline(v=2, col = "blue", lty = 5)
    abline(v=10, col = "red", lty = 5)
dev.off()
}

### check highly conserved ones z > 10
high_cons = plotData$labels[plotData$z> 10]
