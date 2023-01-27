### prep_data.R
## here we process the combined data from TMT consensus cohorts  
## obtained from Eric et al Nat Neuro

if(T){
  rm(list = ls())
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/proteomics/TMT_consensus/Consensus_Banner_ROSMAP/data/")
  
  pacman::p_load(limma,dplyr,tidyr,
                 sva, ## for batch correction
                 imputeLCMD ## package for imputation
                 )

  ## metadata
  datMeta = read.csv("Unified_Consensus_TMT_traits.csv",
                     row.names = 1,skip = 1)
  colnames(datMeta) = gsub("Johnson","",colnames(datMeta))
  colnames(datMeta) = gsub("msex","Sex",colnames(datMeta))
  colnames(datMeta) = gsub("age_death","Age",colnames(datMeta))
  colnames(datMeta) = gsub("braaksc","Braak",colnames(datMeta))
  datMeta$Diagnosis = gsub("Control","CTL",datMeta$Diagnosis)
  datMeta = subset(datMeta,Diagnosis %in%  c("CTL","AsymAD","AD"))
  datMeta$Sex[datMeta$Sex==1] = "male"
  datMeta$Sex[datMeta$Sex==0] = "female"
  datMeta$Age[datMeta$Age == "90+"] = 90
  datMeta$Age = as.numeric(datMeta$Age)
  
  table(datMeta$Diagnosis)
    
  ## load expression data
  datExp = read.csv("2a. Batch and Cohort Corrected Normalized log2 Protein Abundance (0-centered)/2.ConsensusMainArm-BatchAndCohortCorrected_log2_0-centered_Protein_Abundance.csv",
                    row.names = 1)
  
  ## match with expression data 
  idx =intersect(colnames(datExp),rownames(datMeta))
  datExp = datExp[,idx]
  datMeta = datMeta[idx,]
  all(rownames(datMeta)==colnames(exp))
  
  id = which(datMeta$Outlier!="TRUE")
  datMeta = datMeta[id,]
  datExp = datExp[,id]
  
  # datMeta$Cohort[is.na(datMeta$cohort)]=str_split(rownames(datMeta[is.na(datMeta$cohort),]),
  #                                                      pattern = "\\.",simplify = T)[,1]
  datMeta$Batch = factor(datMeta$Cohort)
  datMeta$Diagnosis = factor(datMeta$Diagnosis,levels = c("CTL","AsymAD","AD"))
  
  ## impute missing values ####
  exp = imputeLCMD::impute.QRILC(as.matrix(datExp))
  datExp = exp[[1]]
  
  ## batch correction
  # mod = model.matrix(~ Diagnosis + Sex,datMeta)
  # dt = sva::ComBat(datExp,batch = datMeta$Batch,mod = mod)
  
  ## pca 
  pc = prcomp(t(datExp))
  factoextra::fviz_pca_ind(pc,col.ind = datMeta$Diagnosis,label = "none")
  
  save(file = "C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/data/normalized.Rdata",datExp,datMeta)
  
}
