### prep_data.R
## here we process the combined data from consensus cohorts  
## obtained from Eric et al Nat Med 

if(T){
  rm(list = ls())
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/consensus/")
  
  pacman::p_load(limma,dplyr,tidyr,
                 sva, ## for batch correction
                 imputeLCMD ## package for imputation
                 )

    ## metadata
  datMeta = read.csv("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/proteomics/Erik_consensus/metadata/0.Traits.csv",
                     row.names = 1)
  colnames(datMeta) = gsub("Group","Diagnosis",colnames(datMeta))
  datMeta$Diagnosis = gsub("Control","CTL",datMeta$Diagnosis)
  datMeta = subset(datMeta,Diagnosis %in%  c("CTL","AsymAD","AD"))
  table(datMeta$Diagnosis)
    
  ## load expression data
  datExp = read.csv("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/proteomics/Erik_consensus/raw_data/3.cleanDat.csv",
                    row.names = 1)
  ## match with expression data 
  idx =intersect(colnames(datExp),rownames(datMeta))
  datExp = datExp[,idx]
  datMeta = datMeta[idx,]
  all(rownames(datMeta)==colnames(exp))
  
  datMeta$Batch[is.na(datMeta$Batch)]=str_split(rownames(datMeta[is.na(datMeta$Batch),]),
                                                       pattern = "\\.",simplify = T)[,1]
  datMeta$Batch = factor(datMeta$Batch)
  datMeta$Diagnosis = factor(datMeta$Diagnosis,levels = c("CTL","AsymAD","AD"))
  
  ## impute missing values ####
  exp = imputeLCMD::impute.QRILC(as.matrix(datExp))
  datExp = exp[[1]]
  
  ## batch correction
  mod = model.matrix(~Diagnosis+Sex,datMeta)
  dt = sva::ComBat(datExp,batch = datMeta$Batch,mod = mod)
  
  ## pca 
  pc = prcomp(t(dt))
  factoextra::fviz_pca_ind(pc,col.ind = datMeta$Diagnosis,label = "none")
  
  datMeta$Sex[datMeta$Sex==1] = "male"
  datMeta$Sex[datMeta$Sex==0] = "female"
  save(file = "./data/normalized.Rdata",datExp,datMeta)
  
}
