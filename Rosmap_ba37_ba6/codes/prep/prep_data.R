### prep_data.R
## here we process the combined data from ROSMAP BA6 , BA37  
## obtained from Eric et al Nat Neuro

if(T){
  rm(list = ls())
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/proteomics/TMT_consensus/ROSMAP/ ")
  
pacman::p_load(limma,dplyr,tidyr,
                 sva, ## for batch correction
                 imputeLCMD ## package for imputation
                 )

  ## metadata
  datMeta = as.data.frame(readxl::read_xlsx("0.ROSMAP_BA6+BA37_paired-Traits.xlsx"))
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
  datMeta$Age[is.na(datMeta$Age)]=mean(datMeta$Age,na.rm = T)
  datMeta$Sex[is.na(datMeta$Sex)]=mean(datMeta$Age,na.rm = T)
  rownames(datMeta)= datMeta$batch.channel
  table(datMeta$Diagnosis)
    
  ## load expression data
  datExp = read.csv("2b.BA6+BA37_paired-normalized_relative_abundance-central_tendency_batch_corrected.csv",
                    row.names = 1)
  
  ## match with expression data 
  idx =intersect(colnames(datExp),rownames(datMeta))
  datExp = datExp[,idx]
  datMeta = datMeta[idx,,drop = F]
  rownames(datMeta) = datMeta$batch.channel
  
  ## keep BA37 region (temporal ctx)
  idx = which(datMeta$region=="BA37")
  datMeta = datMeta[idx,]
  datExp = datExp[,idx]
  all(rownames(datMeta)==colnames(datExp))
  
  ## remove NA sex
  idx = which(datMeta$Sex !="NA")
  datMeta = datMeta[idx,]
  datExp = datExp[,idx]
  
  all(rownames(datMeta)==colnames(datExp))
  datMeta$Diagnosis = factor(datMeta$Diagnosis,levels = c("CTL","AsymAD","AD"))
  
  ## impute missing values ####
  exp = imputeLCMD::impute.QRILC(as.matrix(datExp))
  datExp = exp[[1]]
  
  ## batch correction
  # mod = model.matrix(~ Diagnosis + Sex,datMeta)
  # dt = sva::ComBat(datExp,batch = datMeta$Batch,mod = mod)
  
  ## pca 
  pc = prcomp(t(datExp))
  factoextra::fviz_pca_ind(pc,col.ind = datMeta$Diagnosis)
  
  #b21.131N is outlier
  idx = which(rownames(datMeta)=="b21.131N")
  datMeta = datMeta[-idx,]
  datExp = datExp[,-idx]
  pc = prcomp(t(datExp))
  factoextra::fviz_pca_ind(pc,col.ind = datMeta$Diagnosis)
  
  save(file = "C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/Rosmap_ba37_ba6/data/normalized.Rdata",datExp,datMeta)
  
}
