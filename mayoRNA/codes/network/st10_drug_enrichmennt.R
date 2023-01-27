# drug enrichmennt for list of gens from SMR
## here we check enrichment of drugs by rDGIdb v.1.20.0 R package to query DGIdb.

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/analysis_rawCount/disease_specific/AD_PA_CTL/sex_specific/")
library(rDGIdb)
## load smr results
smr = read_delim("../data/gwas/eqtl_results/Brain_cis_eql_smr_results.txt") 

genes = smr$Gene
res = queryDGIdb(genes)
df = res@detailedResults
# df = df[df$FDA==1,]

## save results
write.csv(df,file = "./results/tables/network/smr_drug_targets.csv")
