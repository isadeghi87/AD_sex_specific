## gene enrichment analysis for each module
if(T){
rm(list=ls()); options(stringsAsFactors = F)
library(pacman)
p_load(gprofiler2,ggplot2,dplyr,ggthemes)

setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/replication_proteomics_CSF1/")
attr = read_rds("C:/Users/crgcomu/Desktop/Iman/Brain_meta/annotation/gencode28_attributes.rds")
sex = c("female","male")
res = data.frame()

for (s in sex){
  
  # Loading data ####
  load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
  
  ## filter each module
  for (m in unique(mods)) {
    id = which(mods == m)
    idx = gene[id]
   query = attr$gene_id[match(idx,attr$gene_name)]
   query = query[!is.na(query)]
## Calculate GO enrichment of Disease-associated modules

  go = gost(query,
            organism = "hsapiens",
            significant = FALSE,
            # user_threshold = 0.05,
            correction_method = "fdr",
            domain_scope = "annotated",
            sources = c("KEGG","GO","HPA"))
  
  go = as.data.frame(go$result)
  go = cbind(m,s,go)
  ## order results by p value
  go = go[order(go$p_value),]
  go = go[1:10,]
  res = rbind(res, go)
  
}
}

mod_col = unique(names(mods))
names(mod_col) = unique(mods)

res$term_name = substring(res$term_name, 1, 20)
## plot go for all modules
p = ggplot(res, aes(x = -log(p_value),
                    y = reorder(term_name,log10(p_value)),
           fill = m))+
  geom_bar(stat = "identity",show.legend = F,position = "dodge")+
  facet_wrap(~s,scales = "free_y",ncol=1)+
  scale_fill_manual(values = mod_col)+
 scale_x_log10()+
  labs(x = "-log10(FDR)",y = "")+
  theme_bw(base_size = 9)

p

# ggsave("./results/figures/network/modules_go.pdf",plot = p,
#        width = 15, height = 15)

write.csv(file="./results/tables/network/sex_modules_GO_results.csv", res)


}

