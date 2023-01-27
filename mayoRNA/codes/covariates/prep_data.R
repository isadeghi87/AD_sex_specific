## prepare and normalize the data
#prep_data.R
if(T){
  rm(list=ls())
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/")
  
  ## packages
  pacman::p_load(readr,limma,sva,WGCNA,ggpubr)
  
  #### load data from bulk analysis ####
  #### PA ####
  load("../condition_overlap/disease_specific/PA_metaAnalysis/PA_normalized.Rdata")
  PA_datMeta = datMeta
  rm(datExp)
  
  #### AD ####
  load("../condition_overlap/disease_specific/AD_metaAnalysis/AD_normalized.Rdata")
  AD_datMeta = datMeta
  rm(datExp,datExp.comb,datMeta)
  AD_datMeta$project = "mayo"
  
  #### merge data ####
  datExp = cbind(PA_datExpr,AD_datExpr)
  datmeta = rbind(PA_datMeta,AD_datMeta)

    ### remove frontal & Cerebellum ###
  id = which( datmeta$Study =="Mayo" & datmeta$Brain_Lobe =="Temporal")
  datmeta = datmeta[id,]
  datExp = datExp[,id]
  datmeta$Brain_Lobe = factor(datmeta$Brain_Lobe)
  id = !duplicated(datmeta$Sample)
  
  datmeta = datmeta[id,]
  datExp = datExp[,id]
  
  ## load metadata to add APOE
  meta_apoe = read.csv("../condition_overlap/studies/Mayo/metadata/Mayo_all_Covariates.csv")
  
  datmeta$APOE = meta_apoe$ApoE[match(datmeta$Sample,meta_apoe$Sample)]
  datmeta$APOE = as.factor(datmeta$APOE)
  
  ### keep age>65 & RIN > 7
  id = which(datmeta$Age > 65 & datmeta$RIN > 7)
  datmeta = datmeta[id,]
  datExp = datExp[,id]
  
  ### filter lowly expressed genes
  id = rowSums(cpm(datExp)>0.5)>= 0.3 *ncol(datExp)
  datExp  = datExp[id,]
  
  # dt = sva::ComBat_seq(datExp,batch = datmeta$Brain_Lobe,group = datmeta$Dx)
  ## normalize
  mod <- model.matrix(~0+ Dx + Sex , data = datmeta)
  colnames(mod) = gsub("Dx","",colnames(mod))
  v = voom(datExp, mod,normalize.method = "scale")
  datExp = v$E
  
  ## regress the RIN and Age effect
  mod <- model.matrix(~Dx + Sex + APOE , data = datmeta)
  datExp = removeBatchEffect(datExp,design = mod)
 
   ## Remove outliers based on network connectivity z-scores
  # normadj <- (0.5+0.5*bicor(datExp))^2 ## Calculate connectivity
  # netsummary <- fundamentalNetworkConcepts(normadj)
  # ku <- netsummary$Connectivity
  # z.ku <- (ku-mean(ku))/sqrt(var(ku))
  # outliers = (z.ku < -2)
  # table(outliers)
  # datExp = datExp[,!outliers]
  # datmeta = datmeta[!outliers,]
  
  ## pca to see outliers
  pc = prcomp(t(datExp))
  p = factoextra::fviz_pca_ind(pc,col.ind = datmeta$Dx,
                                label="none"
                                )+
    scale_color_brewer(palette="Set1")+ 
    theme_bw()+
    theme(panel.grid = element_blank())
  
  
  ## keep only protein coding genes ####
  ## read annotation data 
  attr = read_rds("C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/gencode28.rds")
  attr = attr[attr$V3=="gene",]
  attr = attr[match(rownames(datExp),attr$gene_id),]
  attr = subset(attr,gene_type=="protein_coding")
  datExp = datExp[attr$gene_id,]
  
  save(file = "./codes/covariates/normalized_data.Rdata",datExp,datmeta)  
}