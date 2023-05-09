#Simulation comparison MFClust, skm, sgm, R2 cutoff('r' in this script) selected from MFClust

### Simulate data ##
rm(list=ls())
NG = 1000
N = 200

for(l in 1:50){
  set.seed(l)
  Y = matrix(rnorm(N*NG,0,1),nrow = NG,ncol = N)
  
  cl1 = sample(1:3,N,replace = T)
  cl2 = sample(1:3,N,replace = T)
  cl3 = sample(1:3,N,replace = T)
  
  mu.list = c(0.45, 0.47, 0.5, 0.55, 0.6, 0.7, 0.8, 1, 2)
  data.list = list()
  for (i in 1:length(mu.list)) {
    m1 = m2 = m3 = mu.list[[i]]
    
    #View1: 3clusters; 
    Y[1:30,cl1==2] = rnorm(30*sum(cl1==2),m1,1); Y[1:30,cl1==3] = rnorm(30*sum(cl1==3),-m1,1)
    Y[31:60,cl1==1] = rnorm(30*sum(cl1==1),m1,1); Y[31:60,cl1==2] = rnorm(30*sum(cl1==2),-m1,1)
    Y[61:90,cl1==3] = rnorm(30*sum(cl1==3),m1,1); Y[61:90,cl1==1] = rnorm(30*sum(cl1==1),-m1,1)
    
    #View2: 3clusters; 
    Y[91:120,cl2==2] = rnorm(30*sum(cl2==2),m2,1); Y[91:120,cl2==3] = rnorm(30*sum(cl2==3),-m2,1)
    Y[121:150,cl2==1] = rnorm(30*sum(cl2==1),m2,1); Y[121:150,cl2==2] = rnorm(30*sum(cl2==2),-m2,1)
    Y[151:180,cl2==3] = rnorm(30*sum(cl2==3),m2,1); Y[151:180,cl2==1] = rnorm(30*sum(cl2==1),-m2,1)
    
    #View3: 3clusters; 
    Y[181:210,cl3==2] = rnorm(30*sum(cl3==2),m3,1); Y[181:210,cl3==3] = rnorm(30*sum(cl3==3),-m3,1)
    Y[211:240,cl3==1] = rnorm(30*sum(cl3==1),m3,1); Y[211:240,cl3==2] = rnorm(30*sum(cl3==2),-m3,1)
    Y[241:270,cl3==3] = rnorm(30*sum(cl3==3),m3,1); Y[241:270,cl3==1] = rnorm(30*sum(cl3==1),-m3,1)
    
    row.names(Y) = paste0("G",1:NG)
    colnames(Y) = paste0("S",1:N)
    
    data.list[[i]] = t(scale(t(Y)))
  }
  
  clust.lb0 = list(v1 = cl1,
                   v2 = cl2,
                   v3 = cl3)
  gene.lb0 = list(v1 = 1:90, v2=91:180, v3 = 181:270)
  gene.lb0_multi = c(rep(1,90),rep(2,90),rep(3,90),rep(0,NG-270))
  
  save(data.list,clust.lb0,gene.lb0,mu.list,
       file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",l,".RData"))
}

#### MFClust ###################################################################
## r = 0 ##=====================================================================
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
R2.list = c(0)
initial.list = expand.grid(seq(20,60,20),1:10)
colnames(initial.list) = c("kmin","seednum")
mu0.list = c(0.45, 0.47, 0.5, 0.55, 0.6, 0.7, 0.8, 1, 2)

library(parallel)
for(dataIdx in 2:20){
  set.seed(dataIdx)
  print(paste0("Run ",dataIdx,"th data"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  scenario.idx = expand.grid(1:nrow(initial.list), 1:length(R2.list), 1:length(data.list),1:length(V.list))
  full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
  save(full.res,file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.RData"))
  rm(full.res)
}


#Select r_2 cutoff for each mu0 2folds
rm(list=ls())
r2cut.tab.ls = list()
for(dataIdx in 1:20){
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_r20.RData"))
  res.tb.ls = lapply(full.res, "[[","res.tb")
  res.tb = data.frame(do.call(rbind,res.tb.ls))
  res.tbR20.1 = res.tb[res.tb[,"mu0"] >=0.45 & res.tb[,"R2_cutoff"]==0,]
  
  library(dplyr)
  df01 = res.tbR20.1 %>% filter(max_pairARI<0.2) %>%
    group_by(mu0) %>%
    slice(which.max(avgR2_selected_soft_sepV))
  
  r2_cut = c(seq(0.01,0.1,0.01),seq(0.12,0.5,0.02))
  
  r2_cut_thrshld = data.frame(t(sapply(1:nrow(df01), function(i){
    mu0 = df01[i,"mu0"]
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
    ratio.df$mu0 = as.numeric(mu0)
    ratio.df_thrshld2 = ratio.df[which.min(abs(ratio.df$ratio-2)),]
    return(ratio.df_thrshld2)
  })))
  r2_cut_thrshld$mu0 = as.numeric(r2_cut_thrshld$mu0)
  xx = data.frame(mu0 = unique(res.tb$mu0))
  r2cut.tab.ls[[dataIdx]] = left_join(xx,r2_cut_thrshld)
}

avg.tab_r2cut = r2cut.tab.ls[[1]]
for(i in 1:nrow(avg.tab_r2cut)){
  for (j in 1:ncol(avg.tab_r2cut)) {
    avg.tab_r2cut[i,j] = mean(unlist(sapply(r2cut.tab.ls, function(xx) xx[i,j])))
  }
}
avg_r2_cut = round(unlist(avg.tab_r2cut[,"r"]),2)
names(avg_r2_cut) = avg.tab_r2cut[,"mu0"]
save(avg_r2_cut,
     file = "/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_avg_r2_cut_r20_2folds.RData")


### Run for selected r cutoff for each mu0; 2folds===============================
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
source('/home/wez97/MultiViewClust/code/run_func_selectedR2cutoff.R')

load("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_avg_r2_cut_r20_2folds.RData")
V.list = c(3)
R2.list = avg_r2_cut
initial.list = expand.grid(seq(20,60,20),1:10)
colnames(initial.list) = c("kmin","seednum")
mu0.list = c(0.45, 0.47, 0.5, 0.55, 0.6, 0.7, 0.8, 1, 2)

library(parallel)
for(dataIdx in 1:20){
  set.seed(dataIdx)
  print(paste0("Run ",dataIdx,"th data"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  scenario.idx = expand.grid(1:nrow(initial.list), 1:length(data.list),1:length(V.list))
  full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
  save(full.res,file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_select_2folds.RData"))
  rm(full.res)
}

tab.ls = list()
for(dataIdx in 1:20){
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_select_2folds.RData"))
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
    R2cutoff = df$R2_cutoff.0.45[l]
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
     file = "/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_select_2folds_org.RData")


#### SKM ###################################################################
rm(list=ls())
library(Rcpp)
library(S4)
library(parallel)
library(sparcl)
library(caret)

source('/home/wez97/MultiViewClust/code/run_func_skm.R')

load("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_avg_r2_cut_r20_2folds.RData")
R2.list = avg_r2_cut
V = 3

library(parallel)
for(dataIdx in 2:20){
  set.seed(dataIdx)
  print(paste0("Run ",dataIdx,"th data"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  scenario.idx = expand.grid(1:length(R2.list), 1:length(data.list))
  skm.res.ls = mclapply(1:nrow(scenario.idx),run_func,mc.cores = 50)
  save(skm.res.ls,file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/skm/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_select_2folds_skm.RData"))
  rm(skm.res.ls)
}
library(dplyr)
rm(list=ls())
tab.ls = list()
for(dataIdx in 1:20){
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/skm/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_select_2folds_skm.RData"))
  
  
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
      return(c(skm.res.ls[[xx]]$res.tb[1,c("mu0","m")],
               max_pairARI=max_pairARI,avg_pairARI=avg_pairARI))
    }
  })
  res.tb2 = do.call(rbind,res.tb.ls2)
  res.tb2 = res.tb2[!sapply(res.tb2[,1], is.na),]
  res.df = data.frame(res.tb2) %>% group_by(mu0)  %>%
    slice(which.min(max_pairARI))
  
  atab = t(sapply(1:nrow(res.df),function(l){
    mu0 = res.df$mu0[l]
    tb = colMeans(cbind(skm.res.ls[[res.df$m[l]]][["res.tb"]],V = nrow(skm.res.ls[[res.df$m[l]]][["res.tb"]])))
    tb["NG"] = 3*tb["NG"]
    tb
  }))
  df0 = data.frame(mu0 = c(0.45, 0.47, 0.5, 0.55, 0.6, 0.7, 0.8, 1, 2))
  
  yy = merge(df0,atab,"mu0",all = T)
  colnames(yy)[2] = "R2cutoff"
  tab.ls[[dataIdx]] = yy
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
     file = "/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/skm/Simulation_NG1000N200_homoSignal_fullres_select_2folds_skm_org.RData")


#### SGM ###################################################################
rm(list=ls())
library(Rcpp)
library(S4)
library(parallel)
library(sparcl)
library(caret)
library(Brobdingnag)

source('/home/wez97/MultiViewClust/code/run_func_sgm.R')
source("/home/wez97/MultiViewClust/code/sgClust.R")

load("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/mmclust/Simulation_NG1000N200_homoSignal_fullres_avg_r2_cut_r20_2folds.RData")
R2.list = avg_r2_cut
V = 3

library(parallel)
for(dataIdx in 1:20){
  set.seed(dataIdx)
  print(paste0("Run ",dataIdx,"th data"))
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/data/02272023/Simulation_NG1000N200_homoSignal_datalist_",dataIdx,".RData"))
  scenario.idx = expand.grid(1:length(R2.list), 1:length(data.list))
  skm.res.ls = mclapply(1:nrow(scenario.idx),run_func,mc.cores = 50)
  save(skm.res.ls,file = paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/sgm/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_select_2folds_sgm.RData"))
  rm(skm.res.ls)
}

rm(list=ls())
tab.ls = list()
for(dataIdx in 1:20){
  load(paste0("/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/sgm/Simulation_NG1000N200_homoSignal_fullres_",dataIdx,"_select_2folds_sgm.RData"))
  res.tb.ls2 = lapply(1:length(skm.res.ls), function(xx){
    print(xx)
    res = skm.res.ls[[xx]]$sgm.ls
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
      colnames(skm.res.ls[[xx]]$res.tb)[2] = "R2cutoff"
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
     file = "/home/wez97/MultiViewClust/output/Simulation/largeScale/homoSignal/output/02272023/sgm/Simulation_NG1000N200_homoSignal_fullres_select_2folds_sgm_org.RData")

###plot#########################################################################
rm(list=ls())
load("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/output/Simulation_NG1000N200_homoSignal_fullres_select_2folds_org.RData")

library(ggplot2)
p.df.mmclust = data.frame(all.tab[,c("mu0","R2cutoff","v","NG","Cl_ARI","Sensitivity","Specificity")])
p.df.mmclust$model = "MFClust"
avgtb.mmclust = data.frame(avg.tab_R2[,c("mu0","R2cutoff","V","NG","Cl_ARI","Sensitivity","Specificity")])
avgtb.mmclust$model = "MFClust"

library(dplyr)
p.df.mmclust %>% group_by(mu0) %>%
  summarise(xx = median(Cl_ARI))

load("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/output/Simulation_NG1000N200_homoSignal_fullres_select_2folds_skm_org.RData")
p.df.skm = data.frame(all.tab[,c("mu0","R2cutoff","v","NG","Cl_ARI","Sensitivity","Specificity")])
p.df.skm$model = "skm"
avgtb.skm = data.frame(avg.tab_R2[,c("mu0","R2cutoff","V","NG","Cl_ARI","Sensitivity","Specificity")])
avgtb.skm$model = "skm"

load("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/output/Simulation_NG1000N200_homoSignal_fullres_select_2folds_sgm_org.RData")
colnames(all.tab)[2] = colnames(avg.tab_R2)[2] = "R2cutoff"
p.df.sgm = data.frame(all.tab[,c("mu0","R2cutoff","v","NG","Cl_ARI","Sensitivity","Specificity")])
p.df.sgm$model = "sgm"
avgtb.sgm = data.frame(avg.tab_R2[,c("mu0","R2cutoff","V","NG","Cl_ARI","Sensitivity","Specificity")])
avgtb.sgm$model = "sgm"

xx = do.call(rbind,list(avgtb.mmclust, avgtb.skm, avgtb.sgm))
format(xx, digits = 2)

p.df = do.call(rbind,list(p.df.mmclust,p.df.sgm,p.df.skm))
p.df$mu0 = factor(p.df$mu0)
png("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/manuscript/figures/Simulation_select_2folds_clustARI.png",width = 1000,height = 700,res = 200)
ggplot(p.df)+
  geom_boxplot(aes(x = mu0,y = Cl_ARI, fill = model))+
  scale_y_continuous(limits = c(0,1))+
  labs(y = "Cluster ARI",x = bquote(mu))+
  scale_fill_discrete(name = "Model")+
  theme(text = element_text(size=12))
dev.off()

png("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/manuscript/figures/Simulation_select_2folds_Sensitivity.png",width = 900,height = 700,res = 200)
ggplot(p.df)+
  geom_boxplot(aes(x = mu0,y = Sensitivity, fill = model))+
  scale_y_continuous(limits = c(0,1))+
  labs(x = bquote(mu))+
  scale_fill_discrete(name = "Model")+
  theme(text = element_text(size=12))
dev.off()

png("C:/Users/zongw/OneDrive - University of Pittsburgh/Research/MultiViewClust/manuscript/figures/Simulation_select_2folds_Specificity.png",width = 900,height = 700,res = 200)
ggplot(p.df)+
  geom_boxplot(aes(x = mu0,y = Specificity, fill = model))+
  scale_y_continuous(limits = c(0.8,1))+
  labs(x = bquote(mu))+
  scale_fill_discrete(name = "Model")+
  theme(text = element_text(size=12))
dev.off()

