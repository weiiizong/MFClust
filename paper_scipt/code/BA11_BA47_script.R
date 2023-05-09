## BA11_BA47 real data application ### 

#1. V = 2:6, R2cut = 0 =================================================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(tightClust)
library(cluster)
library(MASS)
load("/home/wez97/MultiViewClust/data/BA11_BA47_PreprocessedData.RData")
sourceCpp("/home/wez97/MultiViewClust/code/EM_MFClust_func_C.cpp")
source('/home/wez97/MultiViewClust/code/EM_MFClust_func.R')
source('/home/wez97/MultiViewClust/code/S4_Kmedoids.R')
source('/home/wez97/MultiViewClust/code/run_func_BA11_BA47_r0.R')

std.data = std.data2000
V.list = 2:6
R2.list = c(0)
G.K.list = c(20)
initial.list = expand.grid(c(50,100,500),1:10)
colnames(initial.list) = c("kmin","seednum")
scenario.idx = expand.grid(1:nrow(initial.list), 1:length(R2.list),
                           1:length(G.K.list),1:length(V.list))
library(parallel)
full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
save(full.res,file = "/home/wez97/MultiViewClust/output/RealData/BA11_BA47_03022023_full.res_V23456_r20.RData")

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
     file = "/home/wez97/MultiViewClust/output/RealData/BA11_BA47_03022023_full.res_V23456_avg_r2_cut_r20_2folds.RData") 

#2. V = 2:5, selected R2cut =================================================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(tightClust)
library(cluster)
library(MASS)

load("/home/wez97/MultiViewClust/data/BA11_BA47_PreprocessedData.RData")
sourceCpp("/home/wez97/MultiViewClust/code/EM_MFClust_func_C.cpp")
source('/home/wez97/MultiViewClust/code/EM_MFClust_func.R')
source('/home/wez97/MultiViewClust/code/S4_Kmedoids.R')
source('/home/wez97/MultiViewClust/code/run_func_BA11_BA47_selectedR2cutoff.R')

load("/home/wez97/MultiViewClust/output/RealData/BA11_BA47_03022023_full.res_V23456_avg_r2_cut_r20_2folds.RData")
R2.list = avg_r2_cut
std.data = std.data2000
V.list = 2:5
G.K.list = c(20)
initial.list = expand.grid(c(50,100,500),1:10)
colnames(initial.list) = c("kmin","seednum")
scenario.idx = expand.grid(1:nrow(initial.list), 1:length(G.K.list),1:length(V.list))
library(parallel)
full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
save(full.res,file = "/home/wez97/MultiViewClust/output/RealData/BA11_BA47_03022023_full.res_V2345_selectR2_cutoff_2folds.RData")


###plots=============================================================================
rm(list = ls())
load("~/data/BA11_BA47_PreprocessedData.RData")
std.data = std.data2000
load("~/output/BA11_BA47_03022023_full.res_V2345_selectR2_cutoff_2folds.RData")
res.tb.ls = lapply(full.res, "[[","res.tb")
res.tb = data.frame(do.call(rbind,res.tb.ls))
library(dplyr)
df01 = res.tb %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice(which.max(avgR2_selected_soft_sepV))


p.ls = list()
for (i in 1:nrow(df01)) {
  V = df01[i,"orig_V"]
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
  avg_median_R2 = sapply(R2_V_ls, function(xx) median(xx$sort_r2))
  R2_V_df = do.call(rbind,R2_V_ls)
  R2_V_df$v = factor(R2_V_df$v)
  R2_V_df$genes = factor(R2_V_df$genes,levels = R2_V_df$genes)
  p.ls[[i]] = ggplot(R2_V_df,aes(x = genes, y = sort_r2, color = v))+
    geom_col()+coord_flip()+theme(axis.text.y = element_blank())+
    labs(title = paste0("V=",V))+
    scale_color_manual(values = hcl.colors(V$orig_V,"Zissou1"),
                       labels = paste0(1:V$orig_V,";median_r2=",format(avg_median_R2,digits = 2)))
  
}
library(gridExtra)
grid.arrange(arrangeGrob(grobs= p.ls,ncol=2))

idx = df01$l
library("RColorBrewer")
library(pheatmap)
p.ls = res.tab.ls = list()
for (l in 1:length(idx)) {
  s.idx = idx[l]
  res = full.res[[s.idx]][["res"]]
  pred_GV = res$postGV  
  variableP = full.res[[s.idx]][["metric.tb"]]
  
  
  G.K = full.res[[s.idx]][["res.tb"]]["G.K"]
  V = full.res[[s.idx]][["res.tb"]]["orig_V"]
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
  
  variables = union(row.names(variableP)[apply(variableP, 1, function(x) !all(x>0.05))],c("Age","Brain","Sex"))
  TOD = clinical$TOD
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
  library(tidyr)
  library(ggplot2)
  df = gather(logp.df,variable,log10P,-v)
  p = ggplot(df,aes(x = variable, y = log10P))+
    geom_boxplot()+
    facet_grid(rows = "v",scales = "free_y")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(title = paste0("V=",V))
  p.ls[[l]] = p
  
  variableP = rbind(V,seq(1:res$V),variableP)
  row.names(variableP)[c(1,2)] = c("V","v")
  res.tab.ls[[l]] = variableP
}
library(gridExtra)
grid.arrange(arrangeGrob(grobs= p.ls,ncol=4))
xx = data.frame(t(do.call(cbind,res.tab.ls)))
format(xx,digits = 2)

pairV.ls = list()
for(amu in 1:(length(idx)-1)){
  res1 = full.res[[idx[amu]]][["res"]]
  pred_Glb = apply(res1$postGV,1,function(x) {
    if(all(x == 0)){
      0
    }else{
      which.max(x)
    }
  })
  gene.lb.ls1 = lapply(1:res1$V, function(v){which(pred_Glb == v)})
  
  res2 = full.res[[idx[amu+1]]][["res"]]
  pred_Glb = apply(res2$postGV,1,function(x) {
    if(all(x == 0)){
      0
    }else{
      which.max(x)
    }
  })
  gene.lb.ls2 = lapply(1:res2$V, function(v){which(pred_Glb == v)})
  
  pairV = data.frame(expand.grid(1:length(gene.lb.ls1),1:length(gene.lb.ls2)))
  pairV$weight = sapply(1:nrow(pairV),function(i){
    ridx = pairV[i,1]
    cidx = pairV[i,2]
    
    xx = length(intersect(gene.lb.ls1[[ridx]],gene.lb.ls2[[cidx]]))/min(length(gene.lb.ls1[[ridx]]),length(gene.lb.ls2[[cidx]]))
  })
  pairV$pair = paste0(amu,amu+1)
  pairV.ls[[amu]] = pairV
}
pairV = do.call(rbind,pairV.ls)
pairV.pdf = pairV[pairV$weight != 0,]
V = df01$V
loc = seq(2,12,length.out = length(idx))
nam1 = df01$R2_cutoff.2
vlist = unlist(sapply(V, function(x){
  seq(1:x)
}))
NG = sapply(idx, function(i){
  res = full.res[[i]][["res"]]
  pred_Glb = apply(res$postGV,1,function(x) {
    if(all(x == 0)){
      0
    }else{
      which.max(x)
    }
  })
  gene.lb.ls = sapply(1:res$V, function(v){sum(pred_Glb == v)})
})
pairV.lb = data.frame(x = unlist(sapply(1:length(loc), function(x) rep(loc[x],V[x]))), y = as.numeric(vlist),
                      name = paste0("R2cut=",unlist(sapply(1:length(nam1), function(x) rep(nam1[x],V[x]))),";v",vlist,"\nNG=",unlist(NG)))
library(ggnewscale)
p=ggplot(data = pairV.pdf)
for (i in 1:length(unique(pairV.pdf$pair))) {
  sub.dat = pairV.pdf[pairV.pdf$pair == unique(pairV.pdf$pair)[i],]
  sub.dat$x1 = loc[i]
  sub.dat$x2 = loc[i+1]
  p = p+geom_segment(data = sub.dat,aes(x=x1, y=Var1, xend=x2, yend=Var2,color=weight),size = 2)
}
p = p+scale_color_gradient(low = "white",high = "#0077b6")
p + geom_label(
  data= pairV.lb,
  aes(x = x, y=y),
  label=pairV.lb$name, 
  #nudge_x = 0.25, nudge_y = 0.25, 
  label.size = 1,
  size = 4,
  fill = "#f6bd60")+
  scale_x_continuous(limits = c(0,14))+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  labs(title = paste0("V=",V))

library(mclust)
cl.ls = lapply(1:res$V, function(v){
  apply(res$w_NK_V[[v]],1,which.max)
})

metric.tb = sapply(cl.ls, function(cl){
  Sex = fisher.test(clinical$Sex,cl)$p.value
  Race = fisher.test(clinical$Race,cl)$p.value
  #Cause = fisher.test(clinical$CAUSE,cl)$p.value
  #MANNER = fisher.test(clinical$MANNER,cl)$p.value
  COLSITE = fisher.test(clinical$COLSITE,cl)$p.value
  RAPIDDEATH = fisher.test(clinical$RAPIDDEATH,cl)$p.value
  
  Brain = fisher.test(clinical$Brain,cl)$p.value
  Age = summary(lm(clinical$Age ~ cl))$coefficients[2,4]
  PMI = summary(lm(clinical$PMI ~ cl))$coefficients[2,4]
  pH = summary(lm(clinical$pH ~ cl))$coefficients[2,4]
  RIN = summary(lm(clinical$RIN ~ cl))$coefficients[2,4]
  TOD = summary(lm(clinical$TOD ~ cl))$coefficients[2,4]
  TOD.sin = summary(lm(sin(clinical$TOD) ~ cl))$coefficients[2,4]
  TOD.cos = summary(lm(cos(clinical$TOD) ~ cl))$coefficients[2,4]
  
  return(c(TOD.sin = TOD.sin, TOD.cos = TOD.cos, Sex = Sex, Race = Race,
           COLSITE = COLSITE, RAPIDDEATH = RAPIDDEATH,
           TOD = TOD, Brain = Brain, Age = Age, PMI = PMI, pH = pH, RIN = RIN))
})
colnames(metric.tb) = paste0("view_",1:length(cl.ls))


## V = 5 plot ===========================================================
rm(list = ls())
load("~/data/BA11_BA47_pool_data_03022023.RData")
std.data = std.data2000
load("~/output/BA11_BA47_03022023_full.res_V2345_selectR2_cutoff_2folds.RData")
res.tb.ls = lapply(full.res, "[[","res.tb")
res.tb = data.frame(do.call(rbind,res.tb.ls))
library(dplyr)
df01 = res.tb %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice(which.max(avgR2_selected_soft_sepV))
l = df01$l[4]
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

for (i in 1:length(gene.ls)) {
  write.csv(gene.ls[[i]],paste0("~/data/BA11_BA47_V5_v",i,".csv"))
}


variables = c("Age","Brain","pH","PMI","RIN","Sex")
TOD = clinical$TOD
sinT = sin(TOD)
cosT = cos(TOD)


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
                         "cluster_v5" = factor(cl.ls[[5]]))
  row.names(col_pheno) = colnames(dat)
  annoCol = list("cluster_v1"=c("1"="#D82390", "2"="#E9D2E0"),
                 "cluster_v2"=c("1"="#5E22DF", "2"="#CEC6E1"),
                 "cluster_v3"=c("1"="#1A6DF3", "2"="#C5D8F7"),
                 "cluster_v4"=c("1"="#20D176", "2"="#CCECDC"),
                 "cluster_v5"=c("1"="#F7E526", "2"="#EBE9D6"))
  
  #HC within each cluster to determine order
  cl.ordered = list()
  for(i in 1:length(unique(cl.ls[[v]]))){
    dati = dat[,cl.ls[[v]] == i]
    xx = hclust(dist(t(dati)))
    cl.ordered[[i]] =  colnames(dati)[xx$order]
  }

  dat = apply(dat, 1, standardize)
  
  png(paste0("~/manuscript/figures/BA11_BA47_heatmap_v",v,".png"),width = 800,height = 800,
      res=200)
  pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
               color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
               annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"),annotation_legend = F, annotation_names_col = F,legend = F)
  dev.off()
}
png(paste0("~/manuscript/figures/BA11_BA47_heatmap_legend.png"),width = 1000,height = 1200,
    res=200)
pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
         color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
         annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))
dev.off()

#clinical variable bars
p.ls = list()
for(v in 1:res$V){
  dat = std.data[gene.ls[[v]],]
  col_pheno = data.frame("Cluster" = factor(cl.ls[[v]]),
                         "Age" = clinical$Age,
                         "Brain" = factor(clinical$Brain),
                         "pH" = clinical$pH,
                         "PMI" = clinical$PMI,
                         "RIN" = clinical$RIN,
                         "Sex" = factor(clinical$Sex))
  row.names(col_pheno) = colnames(dat)
  annoCol = list("Cluster"=c("1"="#D82390", "2"="#E9D2E0"),
                 "Brain"=c("B11"="#5E22DF", "B47"="#CEC6E1"),
                 "Sex"=c("F"="#1A6DF3", "M"="#C5D8F7"))
  
  #HC within each cluster to determine order
  cl.ordered = list()
  for(i in 1:length(unique(cl.ls[[v]]))){
    dati = dat[,cl.ls[[v]] == i]
    xx = hclust(dist(t(dati)))
    cl.ordered[[i]] =  colnames(dati)[xx$order]
  }
  
  dat = apply(dat, 1, standardize)
  
  png(paste0("~/manuscript/figures/BA11_BA47_heatmap_v",v,"_topClinical.png"),width = 900,height = 850,
      res=200)
  pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
           color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
           annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"),annotation_legend = T, annotation_names_col = T,legend = F)
  dev.off()
}
png(paste0("~/manuscript/figures/BA11_BA47_heatmap_legend_topClinical.png"),width = 1000,height = 2000,
    res=200)
pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA, 
         color=greenred(75),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
         annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))
dev.off()

metric.tb = lapply(cl.ls, function(cl){
  Age = summary(lm(clinical$Age ~ cl))$coefficients[2,4]
  Brain = fisher.test(clinical$Brain,cl)$p.value
  pH = summary(lm(clinical$pH ~ cl))$coefficients[2,4]
  PMI = summary(lm(clinical$PMI ~ cl))$coefficients[2,4]
  RIN = summary(lm(clinical$RIN ~ cl))$coefficients[2,4]
  Sex = fisher.test(clinical$Sex,cl)$p.value
  sinT = summary(lm(sin(clinical$TOD) ~ cl))$coefficients[2,4]
  cosT = summary(lm(cos(clinical$TOD) ~ cl))$coefficients[2,4]
  
  return(c(Age = Age, Brain = Brain, pH = pH, PMI = PMI, RIN = RIN, Sex = Sex, 
           sinT = sinT, cosT = cosT))
})


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
p.ls = list()
library(ggpubr)

for(v in 1:5){
  metric.tb.v = metric.tb[[v]]
  p.df = logp.df[logp.df$v == v,]
  df = gather(p.df,variable,log10P,-v)
  #df$sigVar_cat = sapply(1:nrow(df), function(xx) metric.tb.v[df$variable[xx]]<0.05)
  df$sigVar = sapply(1:nrow(df), function(xx) -log10(metric.tb.v[df$variable[xx]]))
  df$sigVar[df$sigVar>6] = 6
  ylim = ifelse(v == 4, 190,30)
  p.ls[[v]] = ggplot(df,aes(x = variable, y = log10P,fill = sigVar))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))+
    #scale_fill_continuous(name = bquote(Cluster~-log[10]~P))+
    scale_fill_gradient(
      limits = c(0,6),
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
png("~/manuscript/figures/BA11_BA47_boxplot_V5.png",width = 1200,height = 900,
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
    limits = c(0,6),
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
png("~/manuscript/figures/BA11_BA47_boxplot_V5_legend.png",res=150)
as_ggplot(mylegend)
dev.off()

library(tidyr)
library(ggplot2)
df = gather(logp.df,variable,log10P,-v)
p = ggplot(df,aes(x = variable, y = log10P))+
  geom_boxplot()+
  facet_wrap(~v,scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title = paste0("V=",V))
p


# V=2-5====================================
rm(list = ls())
load("~/data/BA11_BA47_pool_data_03022023.RData")
std.data = std.data2000
load("~/output/BA11_BA47_03022023_full.res_V2345_selectR2_cutoff_2folds.RData")
res.tb.ls = lapply(full.res, "[[","res.tb")
res.tb = data.frame(do.call(rbind,res.tb.ls))
library(dplyr)
df01 = res.tb %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice(which.max(avgR2_selected_soft_sepV))
grid.p = list()
## v=2 plot
l = df01$l[1]
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

variables = c("Age","Brain","pH","PMI","RIN","Sex")
TOD = clinical$TOD
sinT = sin(TOD)
cosT = cos(TOD)

library(mclust)
cl.ls = lapply(1:res$V, function(v){
  apply(res$w_NK_V[[v]],1,which.max)
})

metric.tb = lapply(cl.ls, function(cl){
  Age = summary(lm(clinical$Age ~ cl))$coefficients[2,4]
  Brain = fisher.test(clinical$Brain,cl)$p.value
  pH = summary(lm(clinical$pH ~ cl))$coefficients[2,4]
  PMI = summary(lm(clinical$PMI ~ cl))$coefficients[2,4]
  RIN = summary(lm(clinical$RIN ~ cl))$coefficients[2,4]
  Sex = fisher.test(clinical$Sex,cl)$p.value
  sinT = summary(lm(sin(clinical$TOD) ~ cl))$coefficients[2,4]
  cosT = summary(lm(cos(clinical$TOD) ~ cl))$coefficients[2,4]
  
  return(c(Age = Age, Brain = Brain, pH = pH, PMI = PMI, RIN = RIN, Sex = Sex, 
           sinT = sinT, cosT = cosT))
})


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
p.ls = list()
library(ggpubr)

for(v in 1:res$V){
  metric.tb.v = metric.tb[[v]]
  p.df = logp.df[logp.df$v == v,]
  df = gather(p.df,variable,log10P,-v)
  #df$sigVar_cat = sapply(1:nrow(df), function(xx) metric.tb.v[df$variable[xx]]<0.05)
  df$sigVar = sapply(1:nrow(df), function(xx) -log10(metric.tb.v[df$variable[xx]]))
  df$sigVar[df$sigVar>6] = 6
  ylim = ifelse(v == 4, 190,30)
  p.ls[[v]] = ggplot(df,aes(x = variable, y = log10P,fill = sigVar))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))+
    #scale_fill_continuous(name = bquote(Cluster~-log[10]~P))+
    scale_fill_gradient(
      limits = c(0,6),
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
grid.p[[1]] = grid.arrange(arrangeGrob(grobs= p.ls,nrow = 5))


## v=3 plot
l = df01$l[2]
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

variables = c("Age","Brain","pH","PMI","RIN","Sex")
TOD = clinical$TOD
sinT = sin(TOD)
cosT = cos(TOD)

library(mclust)
cl.ls = lapply(1:res$V, function(v){
  apply(res$w_NK_V[[v]],1,which.max)
})

metric.tb = lapply(cl.ls, function(cl){
  Age = summary(lm(clinical$Age ~ cl))$coefficients[2,4]
  Brain = fisher.test(clinical$Brain,cl)$p.value
  pH = summary(lm(clinical$pH ~ cl))$coefficients[2,4]
  PMI = summary(lm(clinical$PMI ~ cl))$coefficients[2,4]
  RIN = summary(lm(clinical$RIN ~ cl))$coefficients[2,4]
  Sex = fisher.test(clinical$Sex,cl)$p.value
  sinT = summary(lm(sin(clinical$TOD) ~ cl))$coefficients[2,4]
  cosT = summary(lm(cos(clinical$TOD) ~ cl))$coefficients[2,4]
  
  return(c(Age = Age, Brain = Brain, pH = pH, PMI = PMI, RIN = RIN, Sex = Sex, 
           sinT = sinT, cosT = cosT))
})


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
logp.df2 = logp.df
logp.df2$v[logp.df$v == 2] = 3
logp.df2$v[logp.df$v == 3] = 2

p.ls = list()
library(ggpubr)
for(v in 1:res$V){
  metric.tb.v = metric.tb[[v]]
  p.df = logp.df2[logp.df2$v == v,]
  df = gather(p.df,variable,log10P,-v)
  #df$sigVar_cat = sapply(1:nrow(df), function(xx) metric.tb.v[df$variable[xx]]<0.05)
  df$sigVar = sapply(1:nrow(df), function(xx) -log10(metric.tb.v[df$variable[xx]]))
  df$sigVar[df$sigVar>6] = 6
  ylim = ifelse(v == 4, 190,30)
  p.ls[[v]] = ggplot(df,aes(x = variable, y = log10P,fill = sigVar))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))+
    #scale_fill_continuous(name = bquote(Cluster~-log[10]~P))+
    scale_fill_gradient(
      limits = c(0,6),
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
grid.p[[2]] = grid.arrange(arrangeGrob(grobs= p.ls,nrow = 5))

## v=4 plot
l = df01$l[3]
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

variables = c("Age","Brain","pH","PMI","RIN","Sex")
TOD = clinical$TOD
sinT = sin(TOD)
cosT = cos(TOD)

library(mclust)
cl.ls = lapply(1:res$V, function(v){
  apply(res$w_NK_V[[v]],1,which.max)
})

metric.tb = lapply(cl.ls, function(cl){
  Age = summary(lm(clinical$Age ~ cl))$coefficients[2,4]
  Brain = fisher.test(clinical$Brain,cl)$p.value
  pH = summary(lm(clinical$pH ~ cl))$coefficients[2,4]
  PMI = summary(lm(clinical$PMI ~ cl))$coefficients[2,4]
  RIN = summary(lm(clinical$RIN ~ cl))$coefficients[2,4]
  Sex = fisher.test(clinical$Sex,cl)$p.value
  sinT = summary(lm(sin(clinical$TOD) ~ cl))$coefficients[2,4]
  cosT = summary(lm(cos(clinical$TOD) ~ cl))$coefficients[2,4]
  
  return(c(Age = Age, Brain = Brain, pH = pH, PMI = PMI, RIN = RIN, Sex = Sex, 
           sinT = sinT, cosT = cosT))
})


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
logp.df2 = logp.df
logp.df2$v[logp.df$v == 2] = 3
logp.df2$v[logp.df$v == 3] = 2

p.ls = list()
library(ggpubr)
for(v in 1:res$V){
  metric.tb.v = metric.tb[[v]]
  p.df = logp.df2[logp.df2$v == v,]
  df = gather(p.df,variable,log10P,-v)
  #df$sigVar_cat = sapply(1:nrow(df), function(xx) metric.tb.v[df$variable[xx]]<0.05)
  df$sigVar = sapply(1:nrow(df), function(xx) -log10(metric.tb.v[df$variable[xx]]))
  df$sigVar[df$sigVar>6] = 6
  ylim = 30
  p.ls[[v]] = ggplot(df,aes(x = variable, y = log10P,fill = sigVar))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))+
    #scale_fill_continuous(name = bquote(Cluster~-log[10]~P))+
    scale_fill_gradient(
      limits = c(0,6),
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
grid.p[[3]] = grid.arrange(arrangeGrob(grobs= p.ls,nrow = 5))

## v=5 plot
l = df01$l[4]
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

variables = c("Age","Brain","pH","PMI","RIN","Sex")
TOD = clinical$TOD
sinT = sin(TOD)
cosT = cos(TOD)

library(mclust)
cl.ls = lapply(1:res$V, function(v){
  apply(res$w_NK_V[[v]],1,which.max)
})

metric.tb = lapply(cl.ls, function(cl){
  Age = summary(lm(clinical$Age ~ cl))$coefficients[2,4]
  Brain = fisher.test(clinical$Brain,cl)$p.value
  pH = summary(lm(clinical$pH ~ cl))$coefficients[2,4]
  PMI = summary(lm(clinical$PMI ~ cl))$coefficients[2,4]
  RIN = summary(lm(clinical$RIN ~ cl))$coefficients[2,4]
  Sex = fisher.test(clinical$Sex,cl)$p.value
  sinT = summary(lm(sin(clinical$TOD) ~ cl))$coefficients[2,4]
  cosT = summary(lm(cos(clinical$TOD) ~ cl))$coefficients[2,4]
  
  return(c(Age = Age, Brain = Brain, pH = pH, PMI = PMI, RIN = RIN, Sex = Sex, 
           sinT = sinT, cosT = cosT))
})


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
logp.df2 = logp.df

p.ls = list()
library(ggpubr)
for(v in 1:res$V){
  metric.tb.v = metric.tb[[v]]
  p.df = logp.df2[logp.df2$v == v,]
  df = gather(p.df,variable,log10P,-v)
  #df$sigVar_cat = sapply(1:nrow(df), function(xx) metric.tb.v[df$variable[xx]]<0.05)
  df$sigVar = sapply(1:nrow(df), function(xx) -log10(metric.tb.v[df$variable[xx]]))
  df$sigVar[df$sigVar>6] = 6
  ylim = ifelse(v == 4, 190,30)
  p.ls[[v]] = ggplot(df,aes(x = variable, y = log10P,fill = sigVar))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title = paste0("View ",v," (NG=",length(gene.ls[[v]]),")"))+
    #scale_fill_continuous(name = bquote(Cluster~-log[10]~P))+
    scale_fill_gradient(
      limits = c(0,6),
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
grid.p[[4]] = grid.arrange(arrangeGrob(grobs= p.ls,nrow = 5))
png("~/manuscript/figures/BA11_BA47_boxplot_V2345.png",width = 1750,height = 2000,
    res=150)
grid.arrange(grobs = grid.p,ncol = 4)
dev.off()

