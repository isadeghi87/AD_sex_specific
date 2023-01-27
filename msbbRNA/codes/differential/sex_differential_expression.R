### sex_differential_expression.R
## Here we calculate sex-specific differential expression 
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")

library(reshape); 
library(ggplot2);library(biomaRt); library(limma); 
library(readr); library(ggplot2); library(purrr); 
library(stringr); library(dplyr);library(tidyr)
library(statmod);library(RColorBrewer);
library(edgeR); library(ggthemes);library(ggpubr)

#### load data from bulk analysis ####
load("../bulk_analysis/codes/covariates/ad_pa_general.Rdata")

### filter genes
id = rowSums(cpm(datExp)>0.5)>= 0.3 *ncol(datExp)
datExp  = datExp[id,]

## keep only protein coding genes ####
## read annotation data 
attr = read_rds("C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/gencode28.rds")
attr = attr[attr$V3=="gene",]
attr = attr[match(rownames(datExp),attr$gene_id),]
attr = subset(attr,gene_type=="protein_coding")
datExp = datExp[attr$gene_id,]

## normalize
mod <- model.matrix(~ Dx + RIN + Age, data = datmeta)
colnames(mod) = gsub("Dx","",colnames(mod))
v = voom(datExp, mod, normalize.method = "scale")
datExp = v$E
save(file = "./codes/covariates/normalized_data.Rdata",datExp,datmeta,attr)

## check distribution of data 
nsamples <- ncol(datExp)
col <- brewer.pal(nsamples, "Paired")

## samples for each category
p = ggpubr::ggtexttable(table(datmeta$Dx,datmeta$Sex))
ggsave("./results/figures/covariates/bulkRNA_sample_num.pdf",
       p,width = 3,height = 4)

#### differential expression ####
sex = unique(datmeta$Sex)

stats = data.frame()
# stratify for sex
  for (s in sex) {
    idx = which(datmeta$Sex == s)
    Meta = datmeta[idx,]
    Expr = datExp[,idx]
    mod <- model.matrix(~ Dx, data = Meta)
    colnames(mod) = gsub("Dx","",colnames(mod))
    fit <- lmFit(Expr, mod)
    efit <- eBayes(fit, trend = T, robust = T)
    
    ## PA stats
    meta_pa <- topTable(efit, coef = 2, number = Inf, sort.by = "none")
    meta_pa <- meta_pa[,c(1,4,5)]
    meta_pa$sex = s
    meta_pa$Dx = "PA"
    meta_pa$gene_id = rownames(meta_pa)
    rownames(meta_pa) = NULL
    
    ## AD stats
    meta_ad <- topTable(efit, coef = 3, number = Inf, sort.by = "none")
    meta_ad <- meta_ad[,c(1,4,5)]
    meta_ad$sex = s
    meta_ad$Dx = "AD"
    meta_ad$gene_id = rownames(meta_ad)
    rownames(meta_ad) = NULL
    
    stats = rbind(stats,meta_pa,meta_ad)

  }

colnames(stats)[3] = "fdr" 
write.csv(stats, "./results/tables/AD_PA_sex_sumstats.csv")

# spread the data
sig = subset(stats, fdr < 0.05 & abs(logFC)>1)
sig2 = subset(stats, P.Value < 0.05 & abs(logFC)>1)

p1 = ggpubr::ggtexttable(table(sig$Dx,sig$sex))
p2 = ggpubr::ggtexttable(table(sig2$Dx,sig2$sex))

p = ggarrange(p1,p2,labels = c("fdr",
                           "p.value"))
p
ggsave("./results/figures/differential/differential_DEs.pdf",
       p, width = 4,height = 2)

### volcano plot 
p = ggplot(sig2, aes(x = logFC,y = -log10(P.Value)))+
  geom_point(aes(fill = Dx), shape = 21, size = 2)+
  scale_color_manual(values = c("green","red"))+
  facet_grid(Dx ~ sex)+
  theme_bw()+
  theme_base(base_size = 12)

p

ggsave("./results/figures/differential/volcano_plot.pdf",
       width = 10,height = 7)


#### overlap ####
library(venn)
library(GOplot)

venn_plots = list() 
venn_res = list() 

for(i in c(1,2,5,6)){
  j = i+2
  dat1 = deg_list[[i]] %>% dplyr::select(id,logFC)
  dat2 = deg_list[[j]] %>% dplyr::select(id,logFC)
  p = GOVenn(dat1,dat2, label = c("male","female"),plot = F,
             circle.col =c("#185ADB","#FF005C"), 
             lfc.col = c("#4AA96C","#FDCA40","#D62AD0"),# color for directions (up/down)
  )
  n = names(deg_list)[i]
  venn_plots[[n]]=p$plot
  venn_res[[n]] = p$table
  
}

p  = ggarrange(plotlist = venn_plots, 
               ncol = 2, nrow = 2,
               common.legend = T,legend = "bottom",
               labels = c("PA (cerebellum)",'AD (cerebellum)',
                          'PA (temporal)','AD (temporal)'),
               font.label = list(size = 14,face = "plain") 
)
ggsave("./results/overlap/venn_plots.pdf",plot = p,
       width =10 , height = 7)


# 
# ### genes with same direction in both PA and AD in each lobe ####
# 
# gene_dir_list = list()
# 
# for(i in c("up", "down")){
#   # for(s in c("male","female")){
#   for(r in lobes){
#     dat = subset(sig_protein, region ==r&
#                    # sex == s&
#                    direction == i)
#     ids = table(dat$id) %>% sort(decreasing = T)
#     ids = ids[ids ==4]
#     nom = paste(r,i,sep = "_")
#     gene_dir_list[[nom]] = names(ids)
#     # }
#   } 
# }
# 
# 
# ### search for genes up or down in AD but not in PA ####
# dx = "AD"
# 
# ad_unique_deg = list()
# 
# for(dir in c("up","down")){
#   for(r in lobes){
#     for(s in c("male","female")){
#       
#       dat = subset(sig_protein, region ==r & sex == s & direction == dir )
#       
#       ids = table(dat$id) %>% sort(decreasing = T)
#       ids = ids[ids <2]
#       dat = dat[dat$id %in% names(ids) & dat$Dx =="AD",]
#       # print(dat$id)
#       nom = paste(r,s,dx,dir,sep = "_")
#       ad_unique_deg[[nom]] = dat$id
#     }
#   }
# }  
# 
# ## table of unique deg in AD
# ad_table = as.data.frame(matrix(NA, ncol = length(ad_unique_deg),
#                                 nrow = 1))
# colnames(ad_table) = names(ad_unique_deg)
# for(i in colnames(ad_table)){
#   ad_table[,i] = length(ad_unique_deg[[i]])
# }
# 
# #### overlap based on sex ####
# ## same direction in both sexes, stratify by Dx and region
# 
# sex_overlap_dx = list()
# 
# for(i in c("up", "down")){
#   for(r in lobes){
#     for(dx in c("AD","PA")){
#       
#       dat = subset(sig_protein, region ==r&
#                      Dx == dx&
#                      direction == i)
#       ids = table(dat$id) %>% sort(decreasing = T)
#       ids = ids[ids >1]
#       nom = paste(r,dx,i,sep = "_")
#       sex_overlap_dx[[nom]] = names(ids)
#     }
#   } 
# }
# 
# ### table for overlap of both sex 
# sex_table = as.data.frame(matrix(NA, ncol = length(sex_overlap_dx),
#                                  nrow = 1))
# colnames(sex_table) = names(sex_overlap_dx)
# for(i in colnames(sex_table)){
#   sex_table[,i] = length(sex_overlap_dx[[i]])
# }
# 
# 
# mat = sig_protein %>% dplyr::select(id,logFC,region,sex,Dx,gene)  
# mat = pivot_wider(mat, names_from = c(region,Dx,sex),
#                   values_from = logFC)
# mat[is.na(mat)] = 0
# 
# mat2 = pivot_longer(mat,cols = 3:ncol(mat),names_to="var", values_to = "logFC")
# 
# idx = table(sig_protein$id) %>% sort(decreasing = T) %>% head(10) %>% names()
# d = mat2[mat2$id %in% idx,]
# d$gene = factor(d$gene,levels = unique(d$gene))
# d$var = factor(d$var,levels = unique(d$var))
# d$color = ifelse(d$logFC>0, "green","red")
# 
# 
# p = ggplot(d,aes(x = gene, y = var))+
#   geom_tile(aes(fill = logFC))+
#   scale_fill_gradient2( low = "red", high = "darkgreen", na.value="black" )+
#   theme_bw(base_size = 12)+
#   labs(x ="",y= "",title = "Top overlapping genes")+
#   theme(axis.text.x = element_text(angle = 0),
#         plot.title = element_text(hjust = 0.5))
# 
# 
# ggsave("./results/differential_expression/top_overlap_genes.pdf",
#        width = 10, height = 5)

### save data ####
save(file = "./codes/differential/sex_differential_expression.Rdata",sig,sig2,stats)
