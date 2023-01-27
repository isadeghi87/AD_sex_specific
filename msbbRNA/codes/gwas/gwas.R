### gwas.R
setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/")

## libraries

library(ggplot2);library(RRHO);library(tidyverse);
library(corrr);library(corrplot)

plot = T

### load data
load("./codes/differential_expression.Rdata")
sumstats = read.csv("./tables/AD_PA_sumstats.csv",row.names = 1)
sumstats = select(sumstats,contains("logFC"))
colnames(sumstats) = gsub(".logFC","",colnames(sumstats))


getinfo <- c( "ensembl_gene_id",
              "external_gene_name",
              "gene_biotype")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "apr2018.archive.ensembl.org") ## Gencode v28
attr= getBM(
  attributes = getinfo,
  filters = "ensembl_gene_id",
  values = rownames(sumstats),
  mart = mart)

rownames(attr) = attr$ensembl_gene_id


### load GWAS data 

gwas = readxl::read_xlsx(path = "./data/gwas.xlsx",sheet = 18,skip = 5)
gwas = subset(gwas, gwas$`P-VALUE` < 0.000001)

# unique genes
id = unique(gwas$`GWAS gene hit`)
ens_id = attr$ensembl_gene_id[match(id,attr$external_gene_name)]

## subset only gwas data
gwas_exp = sumstats[ens_id,]
gwas_exp = na.omit(gwas_exp)

### overlap of gwas genes with DEGs from each set
library(GeneOverlap)

gwas_overlap = as.data.frame(matrix(ncol = 1, nrow = 8))
rownames(gwas_overlap) = paste(rep(c("Cerebellum","Temporal"),times=1,each=4),
                               rep(c("male","female"),times=1,each=2),
                               rep(c("PA","AD"),times=2,each=1),sep = ".")

## test overlap

for(reg in unique(sig$region)){
  for (dx in unique(sig$condition)) {
    for( sex in unique(sig$sex)){
      ## data for each category
      dat = subset(sig, region == reg & dx == condition & sex == sex)
      idA = unique(dat$id)    
      go.obj = newGeneOverlap(listA = idA,
                              listB = ens_id # gwas IDs
      )
      go.obj = testGeneOverlap(go.obj)
      
      # extract pvalue
      nom = paste(reg,sex,dx,sep = ".")
      gwas_overlap[nom,1] = signif(getPval(go.obj),2)
      
    }
  }
}

gwas_overlap = -log10(gwas_overlap)
gwas_overlap = gwas_overlap[match(colnames(gwas_exp),rownames(gwas_overlap)),]

## annotations
ann = str_split(colnames(gwas_exp),pattern = "\\.",simplify = T)
ann = as.data.frame(ann); colnames(ann) = c("Region","Sex","Condition")
rownames(ann) = colnames(gwas_exp)

## heatmap plot
if(plot){
set.seed(123)

ha = HeatmapAnnotation(Condition = ann$Condition,Region = ann$Region, Sex = ann$Sex)

# p-values for correlation with gwas

ha2 = HeatmapAnnotation( `GWAS overlap \n -log10(p.value)` = gwas_overlap,gp = gpar( fontsize = 10),
                         annotation_name_side = "right",gap = unit(5,"cm"))

pdf("./results/overlap/gwas_heatmap.pdf",width = 6, height = 9)
par(margin(1,0,3,2))
ComplexHeatmap::Heatmap(as.matrix(gwas_exp),name = 'logFC',
                        show_row_names = F,
                        # show_row_dend = F,
                        bottom_annotation= ha,
                        column_names_gp = gpar(fontsize = 8),
                        top_annotation = ha2,
                        height = unit(15,"cm"),
                        width = unit(5,"cm"),column_names_rot = 45)
dev.off()
}


## summary stats
gwas_stats = read.csv("./data/AD_sumstats_Jansenetal_2019sept.txt",sep = "\t")

save(file = "./codes/gwas.Rdata",gwas, gwas_exp,gwas_stats,attr)
