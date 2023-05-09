rm(list=ls())
library(Rcpp)
library(S4)
library(tightClust)
library(cluster)
library(MASS)
load("/home/wez97/MultiViewClust/data/NAc_PreprocessedData.RData")
sourceCpp("/home/wez97/MultiViewClust/code/EM_MFClust_func_C.cpp")
source('/home/wez97/MultiViewClust/code/EM_MFClust_func.R')
source('/home/wez97/MultiViewClust/code/S4_Kmedoids.R')
source('/home/wez97/MultiViewClust/code/run_func_BA11_BA47_r0.R')

std.data = std.data2000
V.list = 2:9
R2.list = c(0)
G.K.list = c(20)
initial.list = expand.grid(c(50,100,500),1:10)
colnames(initial.list) = c("kmin","seednum")
scenario.idx = expand.grid(1:nrow(initial.list), 1:length(R2.list),
                           1:length(G.K.list),1:length(V.list))
library(parallel)
full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
save(full.res,file = "/home/wez97/MultiViewClust/output/RealData/nac_full.res_V23456789_r20.RData")

res.tb.ls = lapply(full.res, "[[","res.tb")
res.tb = data.frame(do.call(rbind,res.tb.ls))
library(dplyr)
df01 = res.tb %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice(which.max(avgR2_selected_soft_sepV))

#Select r_2 cutoff for each mu0
r2_cut = seq(0.01,0.5,0.01)

r2_cut_thrshld = data.frame(t(sapply(1:nrow(df01), function(i){
  orig_V = df01[i,"orig_V"]
  res = full.res[[df01$l[i]]][["res"]]
  pred_Glb = apply(res$postGV,1,function(x) {
    if(all(x == 0)){
      0
    }else{
      which.max(x)
    }
  })
  gene.lb.ls = lapply(1:res$V, function(v){which(pred_Glb == v)})
  
  R2_V_ls = lapply(1:res$V, function(v) {
    r2 = res$R2_V[pred_Glb == v,v]
    df = data.frame(sort_r2 = sort(r2), genes = names(sort(r2)), v = v) 
  })
  ratio.ls = lapply(r2_cut,function(r){
    Vmat = data.frame(t(sapply(1:res$V, function(v){
      if(sum(R2_V_ls[[v]]$sort_r2>r)>=1){
        c(med = median(R2_V_ls[[v]]$sort_r2[R2_V_ls[[v]]$sort_r2>r]), ratio = median(R2_V_ls[[v]]$sort_r2[R2_V_ls[[v]]$sort_r2>r])/r)
      }else{
        c(med = 0, ratio = 0)
      }
    })))
    row.names(Vmat) = c(1:res$V)
    Vmat$r = r
    med = median(Vmat[,"ratio"])
    med.Vmat = Vmat[which.min(abs(Vmat[,"ratio"] - median(Vmat[,"ratio"]))),]
    return(med.Vmat)
  })
  ##median ratio (over V) of median R^2/r is ~2 
  ratio.df = data.frame(do.call(rbind,ratio.ls))
  ratio.df$orig_V = as.numeric(orig_V)
  ratio.df_thrshld2 = ratio.df[which.min(abs(ratio.df$ratio-2)),]
  return(ratio.df_thrshld2)
})))
r2_cut_thrshld$orig_V = as.numeric(r2_cut_thrshld$orig_V)
xx = data.frame(orig_V = unique(df01$orig_V))
avg.tab_r2cut = left_join(xx,r2_cut_thrshld)
avg_r2_cut = round(unlist(avg.tab_r2cut[,"r"]),2)
names(avg_r2_cut) = avg.tab_r2cut[,"orig_V"]
save(avg_r2_cut,
     file = "/home/wez97/MultiViewClust/output/RealData/nac_full.res_V234567_avg_r2_cut_r20_2folds.RData") 

#2. V = 2:6, selected R2cut =================================================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(tightClust)
library(cluster)
library(MASS)

load("/home/wez97/MultiViewClust/data/NAc_PreprocessedData.RData")
sourceCpp("/home/wez97/MultiViewClust/code/EM_MFClust_func_C.cpp")
source('/home/wez97/MultiViewClust/code/EM_MFClust_func.R')
source('/home/wez97/MultiViewClust/code/S4_Kmedoids.R')
source('/home/wez97/MultiViewClust/code/run_func_BA11_BA47_selectedR2cutoff.R')

load("/home/wez97/MultiViewClust/output/RealData/nac_full.res_V234567_avg_r2_cut_r20_2folds.RData")
R2.list = avg_r2_cut
std.data = std.data2000
V.list = 2:7
G.K.list = c(20)
initial.list = expand.grid(c(50,100,500),1:10)
colnames(initial.list) = c("kmin","seednum")
scenario.idx = expand.grid(1:nrow(initial.list), 1:length(G.K.list),1:length(V.list))
library(parallel)
full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
save(full.res,file = "/home/wez97/MultiViewClust/output/RealData/nac_full.res_V234567_selectR2_cutoff_2folds.RData")

## NAc V = 6 plots&tables ===========================================================
rm(list = ls())
load("~/data/NAc_PreprocessedData.RData")

library(biomaRt)
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")
RNAseq = row.names(std.data2000)
output=getBM(attributes=c('ensembl_gene_id','external_gene_name'), 
             filters = 'ensembl_gene_id', 
             values = RNAseq, 
             mart = ensembl)
output.df = data.frame(output)
symbols = output.df$external_gene_name[match(RNAseq,output.df$ensembl_gene_id)]
df = data.frame(RNAseq, symbols)

load("~/output/nac_full.res_V234567_selectR2_cutoff_2folds.RData")
res.tb.ls = lapply(full.res, "[[","res.tb")
res.tb = data.frame(do.call(rbind,res.tb.ls))
library(dplyr)
df01 = res.tb %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice(which.max(avgR2_selected_soft_sepV))
l = df01$l[nrow(df01)]
res = full.res[[l]][["res"]]
pred_GV = res$postGV  
pred_Glb = apply(pred_GV,1,function(x) {
  if(all(x == 0)){
    0
  }else{
    which.max(x)
  }
})
gene.ls = lapply(1:res$V, function(v){
  row.names(std.data2000)[which(pred_Glb == v)]
})

library(openxlsx)
for (i in 1:length(gene.ls)) {
  write.xlsx(data.frame(gene.ls[[i]]),paste0("~/data/NAc_V6_v",i,".xlsx"))
}


symbol.ls = lapply(1:res$V, function(v){
  df$symbols[match(gene.ls[[v]],df$RNAseq)]
})
df5 = data.frame(Ensemble = gene.ls[[5]], Symbol = symbol.ls[[5]])
write.csv(df5, "~/manuscript/figures/snoRNA_in_View5.csv")


rm(list = ls())
load("~/data/NAc_PreprocessedData.RData")
std.data = std.data2000
load("~/output/nac_full.res_V234567_selectR2_cutoff_2folds.RData")
res.tb.ls = lapply(full.res, "[[","res.tb")
res.tb = data.frame(do.call(rbind,res.tb.ls))
library(dplyr)
df01 = res.tb %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice(which.max(avgR2_selected_soft_sepV))
l = df01$l[nrow(df01)]
res = full.res[[l]][["res"]]
pred_GV = res$postGV  
variableP = full.res[[l]][["metric.tb"]]

G.K = full.res[[l]][["res.tb"]]["G.K"]
V = full.res[[l]][["res.tb"]]["orig_V"]
pred_GV = res$postGV  
pred_Glb = apply(pred_GV,1,function(x) {
  if(all(x == 0)){
    0
  }else{
    which.max(x)
  }
})
gene.ls = lapply(1:res$V, function(v){
  row.names(std.data)[which(pred_Glb == v)]
})



library(mclust)
cl.ls = lapply(1:res$V, function(v){
  apply(res$w_NK_V[[v]],1,which.max)
})

library(pheatmap)
library(gplots)
standardize = function(gene.i){
  x = as.numeric(gene.i)
  names(x) = names(gene.i)
  if(sd(x) == 0){
    s.x = (x-mean(x))
  }else{
    s.x = (x-mean(x))/sd(x)
  }
  s.x[which(s.x>2)] = 2
  s.x[which(s.x<(-2))] = -2
  return(s.x)
}
p.ls = list()
for(v in 1:res$V){
  dat = std.data[gene.ls[[v]],]
  col_pheno = data.frame("cluster_v1" = factor(cl.ls[[1]]),
                         "cluster_v2" = factor(cl.ls[[2]]),
                         "cluster_v3" = factor(cl.ls[[3]]),
                         "cluster_v4" = factor(cl.ls[[4]]),
                         "cluster_v5" = factor(cl.ls[[5]]),
                         "cluster_v6" = factor(cl.ls[[6]]))
  row.names(col_pheno) = colnames(dat)
  annoCol = list("cluster_v1"=c("1"="#D82390", "2"="#E9D2E0"),
                 "cluster_v2"=c("1"="#5E22DF", "2"="#CEC6E1"),
                 "cluster_v3"=c("1"="#0861F1", "2"="#D6E0EE"),
                 "cluster_v4"=c("1"="#0EDD34", "2"="#93E6A2", "3"="#D5E6DD"),
                 "cluster_v5"=c("1"="#F7E526", "2"="#EBE9D6"),
                 "cluster_v6"=c("1"="#F35E0D", "2"="#F5DCC5"))
  
  #HC within each cluster to determine order
  cl.ordered = list()
  for(i in 1:length(unique(cl.ls[[v]]))){
    dati = dat[,cl.ls[[v]] == i]
    xx = hclust(dist(t(dati)),method = "ward.D")
    cl.ordered[[i]] =  colnames(dati)[xx$order]
  }
  
  dat = apply(dat, 1, standardize)
  
  png(paste0("~/manuscript/figures/nac_heatmap_v",v,"2.png"),width = 800,height = 800,
      res=200)
  pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
           color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
           annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"),annotation_legend = F, annotation_names_col = F,legend = F)
  dev.off()
}
png(paste0("~/manuscript/figures/nac_heatmap_legend.png"),width = 1000,height = 1200,
    res=200)
pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
         color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
         annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))
dev.off()

p.ls = list()
for(v in 1:res$V){
  dat = std.data[gene.ls[[v]],]
  sex = ifelse(clinical$Sex.Label == "Female","F","M")
  col_pheno = data.frame("Cluster" = factor(cl.ls[[v]]),
                         "pH" = clinical$pH,
                         "RIN" = clinical$RIN,
                         "Sex" = factor(sex))
  row.names(col_pheno) = colnames(dat)
  if(v != 4){
    annoCol = list("Cluster"=c("1"="#D82390", "2"="#E9D2E0"),
                   "Sex"=c("F"="#1A6DF3", "M"="#C5D8F7"))
  }else{
    annoCol = list("Cluster"=c("1"="#D82390","2" = "#DF8EBE", "3"="#E6DCE2"),
                   "Sex"=c("F"="#1A6DF3", "M"="#C5D8F7"))
    
  }

  #HC within each cluster to determine order
  cl.ordered = list()
  for(i in 1:length(unique(cl.ls[[v]]))){
    dati = dat[,cl.ls[[v]] == i]
    xx = hclust(dist(t(dati)))
    cl.ordered[[i]] =  colnames(dati)[xx$order]
  }
  
  dat = apply(dat, 1, standardize)
  
  png(paste0("~/manuscript/figures/nac_heatmap_v",v,"_topClinical.png"),width = 900,height = 850,
      res=200)
  pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
           color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
           annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"),annotation_legend = T, annotation_names_col = T,legend = F)
  dev.off()
}
png(paste0("~/manuscript/figures/nac_heatmap_legend_topClinical2.png"),width = 1000,height = 2000,
    res=200)
pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
         color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
         annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))
dev.off()


variables = c("Age","Sex.Label","PMI","pH","Diagnosis.3Grp","RIN")
TOD = clinical$CorrectedTOD
sinT = sin(TOD)
cosT = cos(TOD)

logp.ls = lapply(1:res$V, function(v){
  sub.std.data = std.data[gene.ls[[v]],,drop=F]
  log.df = data.frame(t(apply(sub.std.data, 1, function(y){
    sapply(variables, function(var){
      -log10(summary(lm(y ~ clinical[,var,drop = T]))$coefficients[2,4])
    })
  })))
  log.df.T = data.frame(t(apply(sub.std.data, 1, function(y){
    sinT.p = -log10(summary(lm(y ~ sinT))$coefficients[2,4])
    cosT.p = -log10(summary(lm(y ~ cosT))$coefficients[2,4])
    return(c(sinT = sinT.p, cosT = cosT.p))
  })))
  log.df = cbind(log.df,log.df.T)
  log.df$v = v
  return(log.df)
})
logp.df = do.call(rbind,logp.ls)
colnames(logp.df) = c("Age","Sex","PMI","pH","Disease","RIN","sinT","cosT","v")

metric.tb = lapply(cl.ls, function(cl){
  Age = summary(lm(clinical$Age ~ cl))$coefficients[2,4]
  Disease = fisher.test(clinical$Diagnosis.3Grp,cl)$p.value
  pH = summary(lm(clinical$pH ~ cl))$coefficients[2,4]
  PMI = summary(lm(clinical$PMI ~ cl))$coefficients[2,4]
  RIN = summary(lm(clinical$RIN ~ cl))$coefficients[2,4]
  Sex = fisher.test(clinical$Sex.Label,cl)$p.value
  sinT = summary(lm(sinT ~ cl))$coefficients[2,4]
  cosT = summary(lm(cosT ~ cl))$coefficients[2,4]
  
  return(c(Age = Age, Disease = Disease, pH = pH, PMI = PMI, RIN = RIN, Sex = Sex, 
           sinT = sinT, cosT = cosT))
})


p.ls = list()
library(ggpubr)

for(v in 1:max(df01$V)){
  metric.tb.v = metric.tb[[v]]
  p.df = logp.df[logp.df$v == v,]
  df = gather(p.df,variable,log10P,-v)
  #df$sigVar_cat = sapply(1:nrow(df), function(xx) metric.tb.v[df$variable[xx]]<0.05)
  df$sigVar = sapply(1:nrow(df), function(xx) -log10(metric.tb.v[df$variable[xx]]))
  df$sigVar[df$sigVar>4] = 4
  if(v == 1){
    ylim = 148
  }else if(v == 3){
    ylim = 30
  }else{
    ylim = 15
  }
  p.ls[[v]] = ggplot(df,aes(x = variable, y = log10P,fill = sigVar))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))+
    #scale_fill_continuous(name = bquote(Cluster~-log[10]~P))+
    scale_fill_gradient(
      limits = c(0,4),
      low = "#C7D1D3",
      high = "#18C3E6",
      space = "Lab",
      na.value = "#C7D1D3",
      guide = "colourbar",
      aesthetics = "fill"
    )+
    #scale_colour_manual(values = c("#a3b18a","#e76f51"))+
    #scale_fill_continuous(limits = c(1.3,6),name = paste0("Cluster -log10P"))+
    scale_y_continuous(limits = c(0,ylim),name = "")+
    scale_x_discrete(name = "")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position = "none",axis.text=element_text(size=11),axis.title=element_text(size=13))
  
}
library(gridExtra)
png("~/manuscript/figures/nac_boxplot_V6.png",width = 1200,height = 900,
    res=150)
grid.arrange(arrangeGrob(grobs= p.ls,ncol=3))
dev.off()
library(ggpubr)
p = ggplot(df,aes(x = variable, y = log10P,fill = sigVar))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  labs(title = paste0("View ",v))+
  #scale_fill_continuous(name = bquote(Cluster~-log[10]~P))+
  scale_fill_gradient(
    limits = c(0,4),
    low = "#C7D1D3",
    high = "#18C3E6",
    space = "Lab",
    na.value = "#C7D1D3",
    guide = "colourbar",
    aesthetics = "fill",
    name = ""
  )+
  scale_y_continuous(limits = c(0,ylim),name = "")+
  scale_x_discrete(name = "")+
  theme(plot.title = element_text(hjust = 0.5))
mylegend = get_legend(p)
png("~/manuscript/figures/nac_boxplot_V6_legend.png",res=150)
as_ggplot(mylegend)
dev.off()
