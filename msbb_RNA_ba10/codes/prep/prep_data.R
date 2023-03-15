### prep_data.R
## RNASEQ obtained MSBB BA10

if(T){
  rm(list = ls())
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/transcriptome/mount_sinai//")
  
  pacman::p_load(limma,dplyr,tidyr,
                 sva, ## for batch correction
                 imputeLCMD ## package for imputation
                 )

  ## metadata
  datMeta = readxl::read_xlsx("0.MSBB.Traits.RNA.xlsx")
  datMeta$Diagnosis = gsub("Control","CTL",datMeta$Diagnosis)
  datMeta = subset(datMeta,Diagnosis %in%  c("CTL","AsymAD","AD"))
  datMeta$Sex[datMeta$Sex==1] = "male"
  datMeta$Sex[datMeta$Sex==0] = "female"
  rownames(datMeta) = datMeta$Specimen.rnaID
  table(datMeta$Diagnosis)
    
  ## load expression data
  datExp = read.csv("3.cleanDat-MSBB-BM10_RNA_all_TRANSPOSED_19070Rx193C.csv",
                    row.names = 1)
  gene = str_split(rownames(datExp),"\\|EN",n = 2,simplify = T)[,1]
  
  ## match with expression data 
  idx = intersect(colnames(datExp),rownames(datMeta))
  datExp = datExp[,idx]
  datMeta = datMeta[idx,]
  all(rownames(datMeta)==colnames(exp))
  
  # datMeta$Cohort[is.na(datMeta$cohort)]=str_split(rownames(datMeta[is.na(datMeta$cohort),]),
  #                                                      pattern = "\\.",simplify = T)[,1]
  datMeta$Diagnosis = factor(datMeta$Diagnosis,levels = c("CTL","AsymAD","AD"))
  
  ## impute missing values ####
  exp = imputeLCMD::impute.QRILC(as.matrix(datExp))
  datExp = exp[[1]]
  
  ## pca 
  pc = prcomp(t(datExp))
  factoextra::fviz_pca_ind(pc,col.ind = datMeta$Diagnosis,label = "none")
  
  save(file = "C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/ROSMAP_RNA/data/normalized.Rdata",datExp,datMeta)
  
}
