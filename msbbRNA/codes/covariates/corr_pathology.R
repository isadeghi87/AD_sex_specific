### correlation with pathological data
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/")

## libraries

library(ggplot2);library(RRHO);library(tidyverse);
library(corrplot);library(ggthemes);library(ggpubr)

### load data
load("./codes/ad_pa_ctl.Rdata")

### filter genes
id = rowSums(cpm(datExp)>0.5)>= 0.3 *ncol(datExp)
datExp  = datExp[id,]

## merge data

all_data = cbind(datmeta,t(datExp))

# data frame to collect results
cor_res = as.data.frame(matrix(NA, ncol = 3,
                               nrow = nrow(datExp)))
colnames(cor_res) = c("Sex","Age","APOE")
rownames(cor_res)= rownames(datExp)
p_res = cor_res

# correlation for each gene and covariate
for(cov in colnames(cor_res)){
  for(id in rownames(datExp)){
    
    mcor =lm(all_data[,id] ~ all_data[,cov], data = all_data )
    cor_res[id,cov] = summary(mcor)$adj.r.squared
    p_res[id, cov] = anova(mcor)$`Pr(>F)`[1]  }
}

### melt data
cor_res2 = pivot_longer(cor_res,cols = 1:3, names_to = "covariate",values_to = "correlation")

## plot
ggplot(cor_res2, aes(x = correlation, y = covariate, fill = covariate))+
  geom_jitter()+
  scale_x_continuous(limits = c(-0.5,0.5))

save(file = "./codes/correlation_pathology.Rdata",all_data,cor_res)
