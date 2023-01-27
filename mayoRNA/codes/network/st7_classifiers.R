## build classifier models for each module
if(T){
  rm(list=ls())
  library(dplyr)
  library(grid)
  library(caret)
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  library(Rtsne)
  library(pROC)
  library(ggrepel)
  library(ggthemes)
  # library(WGCNA);
  library(ComplexHeatmap)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  sex = c("female","male")
  
  res = data.frame()
  topvar = data.frame()
  
  for (s in sex){
   
     # Loading data ####
    load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
    
    ## filter each module
    for (m in unique(mods)) {
      id = which(mods == m)
      datExp = exp[id,]
      
      meta = meta[,"Dx", drop = FALSE]
      data = cbind(meta,t(datExp))
      
      # Partition data into train and test
      set.seed(1234)
      train.idx <- createDataPartition(y = data$Dx, p = 0.5, list = FALSE)[,1]
      X_train <- data[train.idx, ]
      X_test <- data[-train.idx, ]
      
      # Check performance just on the partition (model fitted on train and tested over test)
      # mod_fit <- train(Dx ~. ,  data = X_train, method="rf")
      # confusionMatrix(data = predict(mod_fit, X_test), reference = as.factor(X_test$Dx))
      
      set.seed(42)
      myGrid <- expand.grid(mtry = c(2, 10, 20, 50, 90)) ## Minimal node size; default 1 for classification
      # Perform crossvalidation
      ctrl1 <- trainControl(method = "cv",
                            number = 5,
                            # summaryFunction = twoClassSummary,
                            verboseIter = TRUE,
                            allowParallel = T,
                            savePredictions = TRUE,
                            classProbs = TRUE)
      
      #### random forest model#
      rf_model <- train(Dx ~ ., data = X_train, 
                          method = "rf",
                          tuneGrid=myGrid,
                          trControl = ctrl1)
      
      
      
      ### glmnet model
      glm_model <- train(Dx ~ .,
                         X_train,
                         metric = "ROC",
                         method = "glmnet",
                         tuneGrid = expand.grid(
                           alpha = 0:1,
                           lambda = 0:10/10),
                         trControl = ctrl1)
      
      #### knn model #
      knn_model <- train(Dx ~ .,
                         X_train,
                         metric = "ROC",
                         method = "knn",
                         tuneLength = 20,
                         trControl = ctrl1)
      
      #### svm model #
      svm_model <- train(Dx ~ .,
                         X_train,
                         metric = "ROC",
                         method = "svmRadial",
                         tuneLength = 10,
                         trControl = ctrl1)
      #### naive bayes #
      nb_model <- train(Dx ~ .,
                        X_train,
                        metric = "ROC",
                        method = "naive_bayes",
                        trControl = ctrl1)
     
       ### compare models #
      model_list <- list(glmmet = glm_model,
                         rf = rf_model,
                         knn = knn_model,
                         svm = svm_model,
                         nb = nb_model)
      resamp <- caret::resamples(model_list)
      
      ## choose best model 
      mod.sum = summary(resamp)
      mod.sum = data.frame(accuracy = mod.sum$statistics$Accuracy[,"Mean"],
                       kappa = mod.sum$statistics$Kappa[,"Mean"])
      mod.sum =  mod.sum %>% arrange(desc(accuracy))
      
      top.mod = rownames(mod.sum)[1]
      final_mod = model_list[[grep(top.mod,names(model_list))]]
      pred = predict(final_mod,X_test)
      sens = confusionMatrix(data = pred, 
                             reference = as.factor(X_test$Dx))
      
      # ## choose top important variable
      imp = varImp(final_mod)
      var = as.data.frame(x = imp$importance)
      var = rownames(var)[1:10]
      df = cbind(var,m,s)
      topvar = rbind(topvar,df)
      
      # -- ! CHECK THIS ! --
      accuracy = as.numeric(round(unname(sens$overall["Accuracy"]),2))
      kappa = as.numeric(round(unname(sens$overall["Kappa"]),2))
      p = as.numeric(signif(unname(sens$overall["AccuracyPValue"]),2))
      # sensitivity = round(unname(sens$byClass["Sensitivity"]),2)
      # specificity = round(unname(sens$byClass["Specificity"]),2)
      col = unique(names(mods[mods==m]))
      res = rbind(res,c(m,col,s,accuracy,kappa,p))
    }
  }
  
colnames(res) = c("module","color","sex","accuracy","kappa","p")
  res$p = as.numeric(res$p)
  res$accuracy = as.numeric(res$accuracy)
  res = res[res$module != "M0",]
  write.csv(res,file = "./results/tables/model_accuracy.csv")
  
  topvar = topvar[topvar$module!="M0",]
  topvar$gene_name = attr$gene_name[match(topvar$gene,attr$gene_id)]
  write.csv(topvar,file = "./results/tables/model_top_important_gene.csv")
  
  
  y_lim = c(0,max(-log10(res$p))+1)
  x_lim = c(0,1)
  col = res$color
  names(col) = res$module
  res = res %>% arrange(desc(accuracy))
  
  pp = ggplot(res, aes(x = accuracy, 
                       y = reorder(module,accuracy),
                      fill = module,
                      alpha = p<0.05))+
    geom_point(size = 4,show.legend = F, shape = 21,color = "black")+
    geom_text_repel(aes(label = module),size = 2,show.legend = F,force = T)+
    facet_wrap(~sex,scale= "free",ncol = 2)+
    scale_fill_manual(values = col)+
    labs(y = "")+
    theme_bw(base_size  = 12)+
    theme(strip.background = element_rect(fill = "darkblue"),
          strip.text = element_text(color = "white",face = "bold"),
          axis.text = element_text(size = 8))
  pp
  ggsave("./results/figures/network/classifier_models_sex.pdf",
         width = 8.5,height = 5,plot = pp)
}

