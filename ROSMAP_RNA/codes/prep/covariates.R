# ovariates.R

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/ROSMAP_RNA/")
load("./data/normalized.Rdata")

## plot covariates for each sex ####
sex = unique(datMeta$Sex)
txt = as.data.frame(table(datMeta$Diagnosis,datMeta$Sex))
colnames(txt) = c("dx","sex","freq")
txt$freq[txt$sex=="female"]= -txt$freq

p = ggplot(txt,aes(x = freq,y = dx,fill = sex,group = sex))+
  geom_bar(stat = "identity",
           show.legend = T, 
           width = 0.6,
           color = "black")+
  scale_x_continuous(limits = c(-150,150),
                     breaks = c(-100,0,100),
                     labels = c("100","0","100"))+  
labs(x ="n",
       y="",
       title = "transcriptome subjects")+
  # geom_text(aes(label = Freq),vjust=-0.5)+
  # facet_wrap(~Var2)+
  theme_few()
p
ggsave("./results/figures/covariates/subjects.pdf",p,
       width = 5,height = 3)

for(s in sex){
  iDiagnosis = which(datMeta$Sex == s)
  meta = datMeta[iDiagnosis,]
  exp = datExp[,iDiagnosis]
  
  pdf(paste("./results/figures/covariates/covaritates_",s,".pdf",sep=""), 
      height = 7, 
      width =10)
  par(mfrow = c(3,4))
  
  cols = c("blue","green","red")
  
  # Subjects
  plot(meta$Diagnosis, col=cols, main= "conditions")
  
  # Age
  A= anova(lm(as.numeric(meta$Age) ~ meta$Diagnosis))
  p= A$"Pr(>F)"[1]
  plot(meta$Age ~ meta$Diagnosis, col= cols, 
       main = paste("Age \np=", signif(p,2), sep=""), ylab="year", xlab="")
  
  # Braak
  A= anova(lm(meta$Braak ~ meta$Diagnosis))
  p= A$"Pr(>F)"[1]
  plot(meta$Diagnosis ~ meta$Braak , col= cols, 
       main = paste("Braak \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  # CERAD
  A= anova(lm(meta$CERAD ~ meta$Diagnosis))
  p= A$"Pr(>F)"[1]
  plot(meta$Diagnosis ~ meta$CERAD , col= cols, 
       main = paste("CERAD \np=", signif(p,2), sep=""), ylab="", xlab="")
  
  
  }
  
  dev.off()
}
