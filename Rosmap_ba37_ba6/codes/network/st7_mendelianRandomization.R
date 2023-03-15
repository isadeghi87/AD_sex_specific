## here we compute mendelian randomization using AD GWAS sumstats and
### pQTL from ROSMAP study

rm(list=ls()); options(stringsAsFactors = F)
library(pacman)
p_load(gprofiler2,ggplot2,ggpubr,dplyr,ggthemes,readr)
library(ieugwasr)
library(data.table)
library(MendelianRandomization)
library(TwoSampleMR)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/TMT_consensus/")

## load GWAS and pQTL data
pqtl = read_delim("C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/proteomics/ROSMAP/raw_data/ROSMAP.for_smr.pQTLv2.txt")
gwas = read_delim(file = "C:/Users/crgcomu/Desktop/Iman/Brain_meta/data/gwas/AD_gwas_stats.txt")
colnames(pqtl)[10]= "BETA"

## significant GWAS
sigGwas = gwas[gwas$P<5*10^(-8),]
dim(sigGwas)

sigGwas = gwas
## Next step is to merge pqTL and GWAS
## To this end rename the second column of the CAD data to "rsid" like this.

data.merge = merge(sigGwas, pqtl, by="SNP")
dim(data.merge)

colnames(data.merge) = gsub(".x",".gwas",colnames(data.merge))
colnames(data.merge) = gsub(".y",".pqtl",colnames(data.merge))

## One of the most important steps in pre-processing is to harmonise 
## the effect allele between the two studies. This needs to be done 
## for the effect direction of the betas and the effect allele frequencies.

#**Question E:** Create a TRUE/FALSE vector to indicate if ldl_effectallele is equal to cad_effectallele. How many SNPS have effect alleles that are aligned and how many are not aligned? Tip: Use toupper() or tolower() to change the letters from lower to capital letter or vice-versa.

table(toupper(data.merge$A1.gwas) == data.merge$A1.pqtl)
table(toupper(data.merge$A2.gwas) == data.merge$A2.pqtl)

inconsistent = which(toupper(data.merge$A1.gwas) != data.merge$A1.pqtl & toupper(data.merge$A2.gwas) != data.merge$A2.pqtl)

data.merge[inconsistent,]
data.merge = data.merge[-inconsistent,]


## The harmonization of the effect alleles can be performed as follows:
#
#- IF TRUE: cad_beta = cad_beta_not_aligned and cad_eaf = cad_eaf_not_aligned
#- IF FALSE: cad_beta = -1*cad_beta_not_aligned and cad_eaf = 1-cad_eaf_not_aligned

#Create the new variables cad_beta and cad_eaf and add them to the dataset data.merge. 

# The final step is to prune or clump the SNPs. 
# This step is performed since many MR methods (e.g. IVW) assume that 
# the ratio estimates all provide independent evidence on the causal effect;
# this occurs when the genetic variants are uncorrelated.
# Pruning removes SNPs which are correlated (measured by the squared correlation r2).
# From a group of correlated SNPs it retains the one with the lowest p-value for the exposure.
# Use the function ieugwasr::ld_clump to prune the data. 
# The algorithm needs to know the rs identifier of the genetic variants (labeled as rsid)
# and the $p$-value of the risk factor or exposure (labeled as pval). 
# Rename the following columns accordingly:

colnames(data.merge)[1]="rsid"
colnames(data.merge)[7]="pval"

data.clump = ieugwasr::ld_clump(data.merge)
dim(data.clump)

## Part 2: Mendelian randomization analysis
mr.input = mr_input(bx = data.clump$BETA.gwas, bxse = data.clump$SE.gwas,
                    by = data.clump$BETA.pqtl, byse = data.clump$SE.pqtl,
                    effect_allele = data.clump$A1.gwas,
                    other_allele = data.clump$A2.gwas,eaf = data.clump$EAF,
                    exposure = "GWAS", outcome = "pQTL", snps = data.clump$rsid)

res = mr_ivw(b_exp = data.clump$BETA.gwas,b_out = data.clump$BETA.pqtl,
       se_exp = data.clump$SE.gwas,se_out = data.clump$SE.pqtl 
        )

# **Question B:**  How do you interpret the Heterogeneity test statistic (Cochran's Q)?

# Reply: There is substantial heterogeneity in the IVW MR model, 
# the Q-statistic is > 500 with 77 degrees of freedom, which is 
# much more heterogeneity than expected by chance. This indicates
# potential pleiotropy of the genetic variants used as IVs. 

mr_plot(mr.input, interactive=FALSE)
mr_loo(mr.input)

## Part 3: Robust methods
#perform sensitivity analysis using the function *mr_allmethods()*. Use *mr_mbe()* to apply the mode-based method.

mr_allmethods(mr.input)
mr_plot(mr_allmethods(mr.input,method = "main"))

## Part 4: MR-PRESSO (optional)

# In the final exercise we are going to use the MR-PRESSO package from Verbanck et al 2018 (https://www.nature.com/articles/s41588-018-0099-7) to perform robust MR and outlier detection.
# In order to install the MR-PRESSO package from github (https://github.com/rondolab/MR-PRESSO) use 
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")
library("MRPRESSO")

mrpresso_in =data.frame(ldl_beta,ldl_se,cad_beta,cad_se)

mrpresso_out = mr_presso(BetaOutcome = "cad_beta",BetaExposure = "ldl_beta",SdOutcome = "cad_se",SdExposure = "ldl_se",OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = mrpresso_in, NbDistribution = 1000,  SignifThreshold = 0.05)
```
**Question B:** What is the result from the global test and how do you interpret it?

*Reply: The global test indicates that there is substantial heterogeneity in the data with a residual sum of squares equal to  548.412, which is significantly more than expected under homogenous effects.*
```{r message=FALSE, warning=FALSE}
mrpresso_out$`MR-PRESSO results`$`Global Test`$RSSobs
mrpresso_out$`MR-PRESSO results`$`Global Test`$Pvalue
```
**Question C:** Which genetic variants are detected as outliers?

*Reply: The following indices were identified as outliers: *
```{r message=FALSE, warning=FALSE}
out_indices = mrpresso_out$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
length(out_indices)
rs[out_indices]
```