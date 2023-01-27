#st2_networkParam.R
#---Here we compute the dendrogram and network parameters for the modules
if(T){
  rm(list=ls())
  pacman::p_load(WGCNA,ggpubr,ComplexHeatmap,RColorBrewer,ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/consensus/")
  load("./data/network/parameters/Combined_for_network.Rdata")
  
  enableWGCNAThreads()
  allowWGCNAThreads()
    #load TOM comptued by rWGCNA
    network = list.files("./data/network/parameters/",pattern = paste0("network_signed_"),full.names = T)
    load(network)
    
    geneTree = hclust(1-as.dist(TOM), method="average")
    
    # Iterate WGCNA parameters for robustness -- this takes a while
    colors = vector(mode="list")
    labels = vector(mode="list")
    
    for (minModSize in c(10,15,20)) {
      for (dthresh in c(0.05, 0.1,0.25)) {
        for(ds in c(3:4)) { 
          pam = FALSE
          print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
          
          tree = cutreeHybrid(
            dendro = geneTree,
            minClusterSize = minModSize,
            pamStage = pam, 
            minSplitHeight = dthresh,
            cutHeight = 1,
            deepSplit = ds,
            distM = as.matrix(1-as.dist(TOM)))
          # merged = mergeCloseModules(
          #   exprData = t(voomExp),
          #   colors = tree$labels,
          #   cutHeight=dthresh)
          colors = cbind(colors, labels2colors(tree$labels))
          labels = c(labels, 
                     paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
        }
      }
    }
    
    ## plot different parameters
    pdf("./results/figures/network/WGCNA_diffParams_.pdf", width=18, height=14)
    par(mar= c(5,5,5,5))
    plotDendroAndColors(geneTree, 
                        colors,
                        guideHang = 0.2,
                        groupLabels = labels,
                        addGuide= TRUE,
                        dendroLabels=FALSE,
                        main="Dendrogram",
                        cex.colorLabels=0.9)
    dev.off()
    
    # Finalized Parameters
    power = stringr::str_extract(pattern = "[0-9]++",network)
    power = as.numeric(power)
    wgcna_parameters = list(powers =  power)
    wgcna_parameters$minModSize = 10
    wgcna_parameters$minHeight = 0.1
    wgcna_parameters$bsize = 4000  ##block size needs to be larger than dim(datExpr)[1]
    wgcna_parameters$ds = 4  ##deep split parameter contorls number of modules
    wgcna_parameters$networkType = "signed"    ## using signed networks
    wgcna_parameters$corFnc = "bicor"
    wgcna_parameters$pamStage = FALSE
    
    tree = cutreeHybrid(
      dendro = geneTree,
      minClusterSize = wgcna_parameters$minModSize,
      pamStage = wgcna_parameters$pamStage,
      minSplitHeight = 0.05,
      cutHeight = 1, 
      deepSplit = wgcna_parameters$ds,
      distM = as.matrix(1-as.dist(TOM)))
    
    # merged = mergeCloseModules(exprData= t(datExp),
    # colors = tree$labels,
    # cutHeight = wgcna_parameters$minHeight)
    
    
    # collect colors 
    colors = labels2colors(tree$labels)
    
    # rename modules
    mods = paste0("M",tree$labels)
    n = length(unique(tree$colors))
    names(mods) = colors
    
    # obtain number of genes per module
    colTable = as.data.frame(table(mods))
    
    # colors
    cols = table(names(mods)) %>% sort(decreasing = T) %>% names
    colTable = colTable[order(colTable$Freq,decreasing = T),]
    colTable$color = names(mods)[match(colTable$mods,mods)]
    colTable$mods = factor(colTable$mods,levels = colTable$mods)
    colTable$color = factor(colTable$color,levels = colTable$color)
    
    ## plot module size
    mod_p = ggplot(colTable, aes(x = Freq,
                                 y = mods,
                                 fill = mods))+
      geom_point(
        colour= "black",
        stroke=1,
        size = 4,shape = 21)+
      geom_segment(aes(x= 1, 
                       xend = Freq, 
                       y = mods, 
                       yend = mods),
                   color = "lightgrey",lty=3) +
      scale_fill_manual(values = levels(colTable$color))+
      scale_x_log10()+
      labs(y = "",
           x = "size")+
      theme_few(base_size = 14)+
      theme(legend.position = "none",
            text = element_text(size = 12, color = "black"),
            strip.background = element_rect(fill = "black"),
            strip.text = element_text(color = "white",face = "bold",size=14))
    mod_p
    
    ggsave("./results/figures/network/moduleSize.pdf", plot = mod_p,
           width = 6, height = 4)
    
    ## plot final module
    pdf(file = "./results/figures/network/finalModule.pdf", 
        height = 2, width = 7)
    
    plotDendroAndColors(geneTree,
                        colors,
                        rowTextAlignment = "center", 
                        groupLabels = "Module",
                        main = "",
                        cex.colorLabels = 0.5,
                        addGuide = T,
                        dendroLabels = F)
    
    dev.off()
    
    MEs = moduleEigengenes(expr = t(datExp),
                           tree$labels, 
                           softPower = wgcna_parameters$powers,
                           scale = F)
    kMEtable = signedKME(t(datExp),
                         MEs$eigengenes)
    tableSup = data.frame(kMEtable)
    
    ######check correlation with metadata
    MEs0 = as.data.frame(MEs$eigengenes)
    
    # MEs0 = orderMEs(MEs0)
    MEs0 = MEs0[,colnames(MEs0) != "ME0"]
    colnames(MEs0) = gsub("ME","M",colnames(MEs0))
    
    # merge datMeta and modules eigengene
    ME.datMeta = cbind(datMeta, MEs0)
    
    # a df of correlations
    moduleCors = as.data.frame(matrix(NA, 
                                      nrow = ncol(datMeta), 
                                      ncol = ncol(MEs0)))
    colnames(moduleCors) = colnames(MEs0);
    rownames(moduleCors) = colnames(datMeta)
    
    # covariates
    moduleCors = moduleCors[c("Diagnosis","Sex","Age","CERAD","Braak"),]
    
    #p values
    modulesPval = moduleCors
    
    #calculate linear regression between modules and covariates
    for (mod in colnames(moduleCors)){
      for (var in rownames(moduleCors)){
        mcor =lm(ME.datMeta[,mod] ~ME.datMeta[,var], data = ME.datMeta )
        moduleCors[var,mod] = summary(mcor)$adj.r.squared
        modulesPval[var, mod] = anova(mcor)$`Pr(>F)`[1]
      }
    }
    
    moduleCors = as.matrix(moduleCors)
    modulesPval = as.matrix(modulesPval)
    save(file = paste("./data/network/mod_cov_cor.Rdata",sep=""),moduleCors,modulesPval)
    
    tableSup = cbind(data.frame(mods),gene,kMEtable)
    write.csv(file="./results/tables/network/Modules_kMEtable.csv", tableSup)
    
    save(file = "./data/network/parameters/finalNetwork.RData",
         datExp,
         datMeta,
         gene,
         tableSup,
         geneTree,
         wgcna_parameters,
         colors,
         cols,
         mods,
         MEs,
         kMEtable)
  
  
}
