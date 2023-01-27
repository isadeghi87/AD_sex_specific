### prep_data.R
## here we process the data from ROSMAP-Banner cohorts 
## obtained from Eric et al Nat Neuro 2022

if(T){
  rm(list = ls())
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/rosmap_banner_protein/")
  
  pacman::p_load(reshape,limma,readr,purrr,stringr,dplyr,tidyr,statmod,edgeR,readxl)

    ## metadata
  datMeta = read_excel("./data/data.xlsx",sheet=2,skip=2)
  datMeta = as.data.frame(datMeta)
  colnames(datMeta) = gsub("JohnsonDiagnosis","Diagnosis",colnames(datMeta))
  datMeta = subset(datMeta,!(Diagnosis %in%  c("GIS","Exclude")) & Outlier !="TRUE")
  table(datMeta$Diagnosis)
    
  ## load expression data
  datExp = read_excel("./data/data.xlsx",sheet = 3,skip = 4)
  attr = datExp[,c(2:3)]
  datExp = datExp %>% select(contains(c("uniqueID","rosmap","banner")))
  rownames(datExp) = datExp$UniqueID
  
  ## remove proteins without values
  id = rowSums(is.na(datExp)) != ncol(datExp)-1
  datExp = datExp[id,]
  rownames(datMeta) = datMeta$SpecimenID.internal
  
  ## match with expression data 
  idx =intersect(colnames(datExp),datMeta$SpecimenID.internal)
  datExp = datExp[,idx]
  datMeta = datMeta[idx,]
  all(rownames(datMeta)==colnames(exp))
  
  ## log transformation
  design = model.matrix(~Diagnosis,datMeta)
  
  ## imputation
  
  v = voom(exp,design)
  datExp = v$E
  datMeta$Sex[datMeta$Sex==1] = "male"
  datMeta$Sex[datMeta$Sex==0] = "female"
  save(file = "./codes/tcx_proteomics_normalized.Rdata",datExp,datMeta)
  
  sumstats = data.frame()
  
  for (s in unique(datMeta$Sex)){
    ## differential expression
    mod = model.matrix(~ Diagnosis + AgeAtDeath + PMI ,data = datMeta)
    
    fit <- lmFit(datExp, mod)
    efit <- eBayes(fit, trend = T, robust = T)
    meta <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
    meta <- meta[,c(1,4,5)]
    colnames(meta) = c("logFC", "p.value", "fdr")
    meta$gene = rownames(meta)
    rownames(meta) = NULL
    meta$sex = s
    sumstats = rbind(sumstats,meta)
  }
  write.table(sumstats,"./results/tables/differential/temporal_proteomics_sex_sumstats.txt",col.names = T,row.names = T)
}
