## correlation of modules with covariates
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")
sex = c("female","male")

## load preserved modules
load("./data/network/parameters/finalNetwork.RData")
set = c("M45","M53")
res = data.frame()
eig = MEs$eigengenes
colnames(eig) = gsub("E","",colnames(eig))

set = c("M45","M53")
eig = eig[,set]
datMerged = cbind(datMeta,eig)

cov = c("Braak","CERAD","apoe_genotype")

dat = tidyr::pivot_longer(datMerged,cols = 57:58,names_to = "mod",values_to = "beta")
dat$Braak = as.factor(dat$Braak)

dat$CERAD = as.factor(dat$CERAD)
brak = ggplot(dat,aes(x = Braak, y = beta, fill = Diagnosis))+
  # geom_violin()+
  geom_boxplot(na.rm=T)+
  labs(title = "Braak stages")+
  facet_grid(mod~Sex)+
  theme_bw()
ggsave("./results/figures/network/Braak_mod_sex.pdf",brak,
       width = 6,height = 4)

cerad = ggplot(dat,aes(x = CERAD, y = beta, fill = Diagnosis))+
  geom_boxplot(na.rm = T)+
  labs(title = "CERAD")+
  facet_grid(mod~Sex)+
  theme_bw()
ggsave("./results/figures/network/CERAD_mod_sex.pdf",cerad,
       width = 6,height = 4)

apo = ggplot(dat,aes(x = apoe_genotype, y = beta, fill = Diagnosis))+
  geom_boxplot(na.rm = T)+
  labs(title = "APOE genotype")+
  facet_grid(mod~Sex)+
  theme_bw()
ggsave("./results/figures/network/apoe_mod_sex.pdf",apo,
       width = 6,height = 4)

