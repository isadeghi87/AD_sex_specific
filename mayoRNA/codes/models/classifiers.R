## build classifier models for each module
if(T){
  
  rm(list=ls())
  library(dplyr); library(grid);library(caret)
  library(reshape2);library(ggplot2);library(randomForest);library(magrittr)
  library(Rtsne); library(pROC);library(ggrepel);library(ggthemes)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
  sex = c("female","male")
  load("./codes/covariates/normalized_data.Rdata")
  
  res = data.frame()
  topvar = data.frame()
  
  for (s in sex){
    
    idx = which(datmeta$Sex == s)
    meta = datmeta[idx,]
    exp = datExp[,idx]
    
    ## merge data for training
    dx = meta[,"Dx", drop = FALSE]
    data = cbind(dx,t(exp))
    
    # Partition data into train and test
    set.seed(1234)
    train.idx <- createDataPartition(y = data$Dx, p = 0.5, list = FALSE)[,1]
    X_train <- data[train.idx, ]
    X_test <- data[-train.idx, ]
    
    set.seed(42)
    myGrid <- expand.grid(mtry = c(2, 10, 20, 50, 90)
                          # splitrule = c("gini", "extratrees"),
                          # min.node.size = 1
    ) ## Minimal node size; default 1 for classification
    
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
    var = sortImp(imp,50)
    var = rownames(var)
    df = cbind(var,s)
    topvar = rbind(topvar,df)
    
    # -- ! CHECK THIS ! --
    accuracy = as.numeric(round(unname(sens$overall["Accuracy"]),2))
    kappa = as.numeric(round(unname(sens$overall["Kappa"]),2))
    p = as.numeric(signif(unname(sens$overall["AccuracyPValue"]),2))
    # sensitivity = round(unname(sens$byClass["Sensitivity"]),2)
    # specificity = round(unname(sens$byClass["Specificity"]),2)
    res = rbind(res,cbind(s,accuracy,kappa,p))
    rm(resamp,rf_model,svm_model,glm_model,knn_model)
  }
  
  colnames(res) = c("sex","accuracy","kappa","p")
  res$p = as.numeric(res$p)
  res$accuracy = as.numeric(res$accuracy)
  write.csv(res,file = "./results/tables/conditions_classifier_results.csv")
  topvar$gene_name = attr$gene_name[match(topvar$var,attr$gene_id)]
  write.csv(topvar,file = "./results/tables/conditions_top_important_genes.csv")
  
  gg = ggplot(res, aes(x = sex, 
                       y =accuracy,
                       fill = sex))+
    geom_bar(stat = "identity",width= 0.3,show.legend = F)+
    # scale_fill_manual(values = col)+
    labs(x = "")+
    geom_text(aes(label = paste("p= ",signif(p,2)),y = accuracy+0.05),
              size = 3)+
    theme_bw(base_size  = 12)+
    theme(strip.background = element_rect(fill = "darkblue"),
          strip.text = element_text(color = "white",face = "bold"),
          axis.text = element_text(size = 8))
  gg
  ggsave("./results/figures/models/conditions_classifier_models_sex.pdf",
         width =  5,height = 4,plot = gg)
}

