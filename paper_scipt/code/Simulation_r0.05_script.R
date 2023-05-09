#Supplement comparison MFClust, MFClust_delegate, skm, sgm, r = 0.05

## MFClust r = 0.05 ##=================================================================================
rm(list = ls())
options(stringsAsFactors = F)
library(Rcpp)
library(S4)
library(tightClust)
library(caret)
library(cluster)
library(MASS)
library(parallel)
sourceCpp("/home/wez97/MultiViewClust/code/EM_MFClust_func_C.cpp")
source('/home/wez97/MultiViewClust/code/EM_MFClust_func.R')
source('/home/wez97/MultiViewClust/code/S4_Kmedoids.R')
source('/home/wez97/MultiViewClust/code/run_func_r0.R')

V.list = c(3)
R2.list = c(0.05)
initial.list = expand.grid(seq(20,60,20),1:10)
colnames(initial.list) = c("kmin","seednum")
mu0.list = c(0.45, 0.47, 0.5, 0.55, 0.6, 0.7, 0.8, 1, 2)

library(parallel)
for(dataIdx in 1:20){
  set.seed(dataIdx)
  print(paste0("Run ",dataIdx,"th data"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  scenario.idx = expand.grid(1:nrow(initial.list), 1:length(R2.list), 1:length(data.list),1:length(V.list))
  full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
  save(full.res,file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/mmclust/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.05_mmclust.RData"))
  rm(full.res)
}
rm(list=ls())
tab.ls = list()
for(dataIdx in 1:20){
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/mmclust/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.05_mmclust.RData"))
  res.tb.ls = lapply(full.res, "[[","res.tb")
  res.tb = data.frame(do.call(rbind,res.tb.ls))
  library(dplyr)
  df = res.tb %>% filter(max_pairARI<0.2) %>%
    group_by(mu0) %>%
    slice(which.max(avgR2_selected_soft_sepV) )
  
  idx = df$l
  tab.ls[[dataIdx]] = t(sapply(1:nrow(df),function(l){
    res = full.res[[df$l[l]]][["res"]]
    mu0 = df$mu0[l]
    R2cutoff = df$R2_cutoff[l]
    tb = colMeans(cbind(mu0,R2cutoff,V = nrow(full.res[[df$l[l]]][["metric.tb"]]),full.res[[df$l[l]]][["metric.tb"]]))
    tb["NG"] = 3*tb["NG"]
    tb
  }))
}
avg.tab_R2 = Reduce(`+`, tab.ls)/length(tab.ls)
sd.tab_R2 = avg.tab_R2
for(i in 1:nrow(sd.tab_R2)){
  for (j in 1:ncol(sd.tab_R2)) {
    sd.tab_R2[i,j] = sd(sapply(tab.ls, function(xx) xx[i,j]))
  }
}
all.tab = do.call(rbind, tab.ls)
save(avg.tab_R2,sd.tab_R2,all.tab,
     file = "/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/mmclust/Simulation_NG1000N200_homoSignal_fullres_r20.05_mmclust_org.RData")

## MFClust delegate r = 0.05 ##=================================================================================
rm(list = ls())
options(stringsAsFactors = F)
library(Rcpp)
library(S4)
library(tightClust)
library(caret)
library(cluster)
library(MASS)
library(parallel)
sourceCpp("/home/wez97/MultiViewClust/code/EM_MFClust_delegate_func_C.cpp")
source('/home/wez97/MultiViewClust/code/EM_MFClust_delegate_func.R')
source('/home/wez97/MultiViewClust/code/S4_Kmedoids.R')
source('/home/wez97/MultiViewClust/code/run_func_r0.R')

V.list = c(3)
R2.list = c(0.05)
initial.list = expand.grid(seq(20,60,20),1:10)
colnames(initial.list) = c("kmin","seednum")
mu0.list = c(0.45, 0.47, 0.5, 0.55, 0.6, 0.7, 0.8, 1, 2)

library(parallel)
for(dataIdx in 1:20){
  set.seed(dataIdx)
  print(paste0("Run ",dataIdx,"th data"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  scenario.idx = expand.grid(1:nrow(initial.list), 1:length(R2.list), 1:length(data.list),1:length(V.list))
  full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
  save(full.res,file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/mmclust_delegate/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.05_mmclust_delegate.RData"))
  rm(full.res)
}
rm(list=ls())
tab.ls = list()
for(dataIdx in 1:20){
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/mmclust_delegate/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.05_mmclust_delegate.RData"))
  res.tb.ls = lapply(full.res, "[[","res.tb")
  res.tb = data.frame(do.call(rbind,res.tb.ls))
  library(dplyr)
  df = res.tb %>% filter(max_pairARI<0.2) %>%
    group_by(mu0) %>%
    slice(which.max(avgR2_selected_soft_sepV) )
  
  idx = df$l
  tab.ls[[dataIdx]] = t(sapply(1:nrow(df),function(l){
    res = full.res[[df$l[l]]][["res"]]
    mu0 = df$mu0[l]
    R2cutoff = df$R2_cutoff[l]
    tb = colMeans(cbind(mu0,R2cutoff,V = nrow(full.res[[df$l[l]]][["metric.tb"]]),full.res[[df$l[l]]][["metric.tb"]]))
    tb["NG"] = 3*tb["NG"]
    tb
  }))
}
avg.tab_R2 = Reduce(`+`, tab.ls)/length(tab.ls)
sd.tab_R2 = avg.tab_R2
for(i in 1:nrow(sd.tab_R2)){
  for (j in 1:ncol(sd.tab_R2)) {
    sd.tab_R2[i,j] = sd(sapply(tab.ls, function(xx) xx[i,j]))
  }
}
all.tab = do.call(rbind, tab.ls)
save(avg.tab_R2,sd.tab_R2,all.tab,
     file = "/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/mmclust_delegate/Simulation_NG1000N200_homoSignal_fullres_r20.05_mmclust_delegate_org.RData")

## skm r = 0.05 ##=================================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(parallel)
library(sparcl)
library(caret)

source('/home/wez97/MultiViewClust/code/run_func_skm.R')

V = 3
R2.list = c(0.05)

library(parallel)
for(dataIdx in 2:20){
  set.seed(dataIdx)
  print(paste0("Run ",dataIdx,"th data"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  scenario.idx = expand.grid(1:length(R2.list), 1:length(data.list))
  skm.res.ls = mclapply(1:nrow(scenario.idx),run_func,mc.cores = 50)
  save(skm.res.ls,file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/skm/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.05_skm.RData"))
  rm(skm.res.ls)
}

rm(list=ls())
tab.ls = list()
for(dataIdx in 1:20){
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/skm/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.05_skm.RData"))

  
  res.tb.ls2 = lapply(1:length(skm.res.ls), function(xx){
    print(xx)
    res = skm.res.ls[[xx]]$skm.ls
    if(length(res) < 2){
      return(NA)
    }else{
      Vpairs = combn(length(res),2)
      sample_Clust_ARI = sapply(1:ncol(Vpairs), function(z) {
        mclust::adjustedRandIndex(factor(res[[Vpairs[1,z]]]$Cs),factor(res[[Vpairs[2,z]]]$Cs))
      })
      max_pairARI = max(sample_Clust_ARI)
      avg_pairARI = mean(sample_Clust_ARI)
      return(c(skm.res.ls[[xx]]$res.tb[1,c("mu0","R2cutoff","m")],
               max_pairARI=max_pairARI,avg_pairARI=avg_pairARI))
    }
  })
  res.tb2 = do.call(rbind,res.tb.ls2)
  res.tb2 = res.tb2[!sapply(res.tb2[,1], is.na),]
  res.df = data.frame(res.tb2) %>% group_by(mu0)  %>%
    slice(which.min(max_pairARI))
  
  atab = t(sapply(1:nrow(res.df),function(l){
    mu0 = res.df$mu0[l]
    R2cutoff = res.df$R2cutoff[l]
    tb = colMeans(cbind(skm.res.ls[[res.df$m[l]]][["res.tb"]],V = nrow(skm.res.ls[[res.df$m[l]]][["res.tb"]])))
    tb["NG"] = 3*tb["NG"]
    tb
  }))
  df0 = data.frame(mu0 = c(0.45, 0.47, 0.5, 0.55, 0.6, 0.7, 0.8, 1, 2))
  tab.ls[[dataIdx]] = merge(df0,atab,"mu0",all = T)
}
avg.tab_R2 = sd.tab_R2 = tab.ls[[1]]
for(i in 1:nrow(sd.tab_R2)){
  for (j in 1:ncol(sd.tab_R2)) {
    sd.tab_R2[i,j] = sd(sapply(tab.ls, function(xx) xx[i,j]),na.rm=T)
    avg.tab_R2[i,j] = mean(sapply(tab.ls, function(xx) xx[i,j]),na.rm=T)
    
  }
}
all.tab = do.call(rbind, tab.ls)
save(avg.tab_R2,sd.tab_R2,all.tab,
     file = "/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/skm/Simulation_NG1000N200_homoSignal_fullres_r20.05_skm_org.RData")

## sgm r = 0.05 ##=================================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(parallel)
library(sparcl)
library(caret)
library(Brobdingnag)

source('/home/wez97/MultiViewClust/code/run_func_sgm.R')
source("/home/wez97/MultiViewClust/code/sgClust.R")

V = 3
R2.list = c(0.05)

library(parallel)
for(dataIdx in 2:20){
  set.seed(dataIdx)
  print(paste0("Run ",dataIdx,"th data"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  scenario.idx = expand.grid(1:length(R2.list), 1:length(data.list))
  sgm.res.ls = mclapply(1:nrow(scenario.idx),run_func,mc.cores = 50)
  save(sgm.res.ls,file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/sgm/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.05_sgm.RData"))
  rm(sgm.res.ls)
}

rm(list=ls())
tab.ls = list()
for(dataIdx in 1:20){
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/sgm/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.05_sgm.RData"))
  res.tb.ls2 = lapply(1:length(sgm.res.ls), function(xx){
    print(xx)
    res = sgm.res.ls[[xx]]$sgm.ls
    if(length(res) < 2){
      return(NA)
    }else{
      Vpairs = combn(length(res),2)
      sample_Clust_ARI = sapply(1:ncol(Vpairs), function(i) {
        cl1 = apply(res[[Vpairs[1,i]]]$result$z,1,which.max)
        cl2 = apply(res[[Vpairs[2,i]]]$result$z,1,which.max)
        
        mclust::adjustedRandIndex(factor(cl1),factor(cl2))
      })
      max_pairARI = max(sample_Clust_ARI)
      avg_pairARI = mean(sample_Clust_ARI)
      return(c(sgm.res.ls[[xx]]$res.tb[1,c("mu0","R2cutoff","m")],
               max_pairARI=max_pairARI,avg_pairARI=avg_pairARI))
    }
  })
  res.tb2 = do.call(rbind,res.tb.ls2)
  res.tb2 = res.tb2[!sapply(res.tb2[,1], is.na),]
  res.df = data.frame(res.tb2) %>% group_by(mu0)  %>%
    slice(which.min(max_pairARI))
  
  atab = t(sapply(1:nrow(res.df),function(l){
    mu0 = res.df$mu0[l]
    R2cutoff = res.df$R2cutoff[l]
    tb = colMeans(cbind(sgm.res.ls[[res.df$m[l]]][["res.tb"]],V = nrow(sgm.res.ls[[res.df$m[l]]][["res.tb"]])))
    tb["NG"] = 3*tb["NG"]
    tb
  }))
  df0 = data.frame(mu0 = c(0.45, 0.47, 0.5, 0.55, 0.6, 0.7, 0.8, 1, 2))
  tab.ls[[dataIdx]] = merge(df0,atab,"mu0",all = T)
}
avg.tab_R2 = sd.tab_R2 = tab.ls[[1]]
for(i in 1:nrow(sd.tab_R2)){
  for (j in 1:ncol(sd.tab_R2)) {
    sd.tab_R2[i,j] = sd(sapply(tab.ls, function(xx) xx[i,j]),na.rm=T)
    avg.tab_R2[i,j] = mean(sapply(tab.ls, function(xx) xx[i,j]),na.rm=T)
    
  }
}
all.tab = do.call(rbind, tab.ls)
save(avg.tab_R2,sd.tab_R2,all.tab,
     file = "/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023_r0.05/sgm/Simulation_NG1000N200_homoSignal_fullres_r20.05_sgm_org.RData")

rm(list=ls())
load("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/output/Simulation_NG1000N200_homoSignal_fullres_r20.05_mmclust_org.RData")

library(ggplot2)
p.df.mmclust = data.frame(all.tab[,c("mu0","R2cutoff","v","NG","Cl_ARI","Sensitivity","Specificity")])
p.df.mmclust$model = "MFClust"
avgtb.mmclust = data.frame(avg.tab_R2[,c("mu0","R2cutoff","V","NG","Cl_ARI","Sensitivity","Specificity")])
avgtb.mmclust$model = "MFClust"

load("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/output/Simulation_NG1000N200_homoSignal_fullres_r20.05_mmclust_delegate_org.RData")
p.df.mmclust_delegate = data.frame(all.tab[,c("mu0","R2cutoff","v","NG","Cl_ARI","Sensitivity","Specificity")])
p.df.mmclust_delegate$model = "MFClust_delegate"
avgtb.mmclust_delegate = data.frame(avg.tab_R2[,c("mu0","R2cutoff","V","NG","Cl_ARI","Sensitivity","Specificity")])
avgtb.mmclust_delegate$model = "MFClust_delegate"

load("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/output/Simulation_NG1000N200_homoSignal_fullres_r20.05_skm_org.RData")
p.df.skm = data.frame(all.tab[,c("mu0","R2cutoff","v","NG","Cl_ARI","Sensitivity","Specificity")])
p.df.skm$model = "skm"
avgtb.skm = data.frame(avg.tab_R2[,c("mu0","R2cutoff","V","NG","Cl_ARI","Sensitivity","Specificity")])
avgtb.skm$model = "skm"

load("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/output/Simulation_NG1000N200_homoSignal_fullres_r20.05_sgm_org.RData")
p.df.sgm = data.frame(all.tab[,c("mu0","R2cutoff","v","NG","Cl_ARI","Sensitivity","Specificity")])
p.df.sgm$model = "sgm"
avgtb.sgm = data.frame(avg.tab_R2[,c("mu0","R2cutoff","V","NG","Cl_ARI","Sensitivity","Specificity")])
avgtb.sgm$model = "sgm"

xx = do.call(rbind,list(avgtb.mmclust, avgtb.mmclust_delegate, avgtb.skm, avgtb.sgm))
format(xx, digits = 2)

p.df = do.call(rbind,list(p.df.mmclust,p.df.mmclust_delegate,p.df.sgm,p.df.skm))
p.df$mu0 = factor(p.df$mu0)
png("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/manuscript/figures/Simulation_r0.05_clustARI.png",width = 1100,height = 700,res = 200)
ggplot(p.df)+
  geom_boxplot(aes(x = mu0,y = Cl_ARI, fill = model))+
  scale_y_continuous(limits = c(0,1))+
  labs(y = "Cluster ARI",x = bquote(mu))+
  scale_fill_discrete(name = "Model")+
  theme(text = element_text(size=12))
dev.off()

png("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/manuscript/figures/Simulation_r0.05_Sensitivity.png",width = 1000,height = 700,res = 200)
ggplot(p.df)+
  geom_boxplot(aes(x = mu0,y = Sensitivity, fill = model))+
  scale_y_continuous(limits = c(0,1))+
  labs(x = bquote(mu))+
  scale_fill_discrete(name = "Model")+
  theme(text = element_text(size=12))
dev.off()

png("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/manuscript/figures/Simulation_r0.05_Specificity.png",width = 1000,height = 700,res = 200)
ggplot(p.df)+
  geom_boxplot(aes(x = mu0,y = Specificity, fill = model))+
  scale_y_continuous(limits = c(0.8,1))+
  labs(x = bquote(mu))+
  scale_fill_discrete(name = "Model")+
  theme(text = element_text(size=12))
dev.off()

