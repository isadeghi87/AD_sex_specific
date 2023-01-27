if(T){
  rm(list=ls()); options(stringsAsFactors = F)
  pacman::p_load(WGCNA, ggplot2,reshape, magrittr,ggthemes,
                 ggstatsplot)
  
  setwd("C:/Users/crgcomu/Desktop/Iman/Brain_meta/projects/AD_sex_specific/")
  
  ## rad DE results for modules
  panel = read.csv("./results/tables/network/modules_DE.csv",row.names = 1)
  
  sex = c("female","male")
  csf = c("CSF1","CSF2")
  
  res = stat = data.frame()
  
  for(dt in csf){
    ## CSF data
    load(paste0("./replication_proteomics_",dt,"/codes/",dt,"_proteomics_normalized.Rdata"))
    
    for (s in sex){
      
      id = which(datMeta$Sex==s)
      submeta = datMeta[id,]
      subexp = datExp[,id]
      
      subpanel = panel[panel$sex ==s,]
      #load final network
      load(paste0("./codes/network/parameters/finalNetwork_",s,".RData",sep=""))
      
      for(cell in unique(subpanel$panel)){
        m = subpanel$module[subpanel$panel == cell]
        m = unique(m)
        
        for(mod in unique(m)){
      
      ## intersect the matched genes
      id = attr$gene_name[mods== mod]
      index = which(gene %in%  id)
      if(length(index)>2){
      gene.index = gene[index]
      mat = subexp[index,]
      
      ## add metadata
      mat = as.data.frame(t(mat))
      mat  = cbind(mat,submeta)
      df = pivot_longer(mat,cols = 1:length(index),names_to = "gene",values_to = "exp")
      df = df[,c("gene","exp","Diagnosis")]
      df$dataset = dt
      df$panel = cell
      df$module = mod
      df$sex = s
      res = rbind(res,df)
      
      ## test between groups
      model = lm(exp~Diagnosis, data = df)
      mod.sum = summary(model)
      mod.sum = mod.sum$coefficients
      r = mod.sum[2,1]
      se = mod.sum[2,2]
      pval = mod.sum[2,4]
      stat = rbind(stat,cbind(pval,r,se,dt,mod,cell,s))
      }
        }
    }
    
    
  }
  }
  
  dat = res[res$dataset=="CSF1",]
  ## plot 
  ggplot(dat, aes(x = Diagnosis, y = exp, fill = panel))+
    geom_boxplot(show.legend = F,outlier.shape = NA)+
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    facet_grid(sex~panel,scales = "free",space="free")+
    scale_y_continuous(limits = c(-2,2))+
    scale_fill_manual(values = c(2:7))+
    labs(y = "beta", x = "", title = "") +
    # scale_x_discrete()+
    theme_bw(base_size = 14)+
    theme(legend.position = "none",
          strip.text = element_text(size = 12,face = "bold",color = "white"),
          axis.title = element_text(size = 12,angle = 90),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.background.y = element_rect(fill = "darkblue"),
          strip.background.x = element_rect(fill = "darkgreen"),
    )
  
}

stat$r = as.numeric(stat$r)
stat$pval = as.numeric(stat$pval)
stat$se = as.numeric(stat$se)
df = stat[stat$s==s,]

ggplot(df,aes(x = mod, y = r, fill = cell))+
  geom_bar(stat= "identity",width = 0.7,show.legend = F)+
  facet_grid(dt~cell,scale="free",space="free")+
  geom_text(aes(label = signif(pval,1)),size=2)+ 
  theme_bw()
