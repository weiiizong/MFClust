### Simulate data ##
rm(list=ls())
NG = 1000
N = 200

mu.list = c(0.8,1,1.2)
for (mu in mu.list) {
  data.list = clust.lb0 = list()
  
  for(l in 1:20){
    set.seed(l)
    Y = matrix(rnorm(N*NG,0,1),nrow = NG,ncol = N)
    
    cl1 = sample(1:3,N,replace = T)
    
    Y[1:270,cl1==2] = rnorm(270*sum(cl1==2),mu,1); Y[1:270,cl1==3] = rnorm(270*sum(cl1==3),-mu,1)
    
    row.names(Y) = paste0("G",1:NG)
    colnames(Y) = paste0("S",1:N)
    
    data.list[[l]] = t(scale(t(Y)))
    
    
    clust.lb0[[l]] = list(v1 = cl1)
    gene.lb0 = list(v1 = 1:270)
  }
  save(data.list,clust.lb0,gene.lb0,
       file = paste0("~/data/Simulation_NG1000N200_homoSignal_datalist_",mu,".RData"))
  
}
rm(list = ls())
options(stringsAsFactors = F)
library(Rcpp)
library(S4)
library(tightClust)
library(caret)
library(cluster)
library(MASS)
library(parallel)
sourceCpp("~/code/EM_MFClust_func_C.cpp")
source('~/code/EM_MFClust_func.R')
source('~/code/S4_Kmedoids.R')
source('~/code/run_funcs_fixR2.R')
V.list = c(1:2)
R2.list = c(0)
initial.list = expand.grid(seq(20,60,20),1:5)
colnames(initial.list) = c("kmin","seednum")
for (mu0 in c(0.8,1,1.2)) {
  load(paste0("~/data/Simulation_NG1000N200_homoSignal_datalist_",mu0,".RData"))
  set.seed(12345)
  scenario.idx = expand.grid(1:nrow(initial.list), 1:length(R2.list), 1:length(data.list),1:length(V.list))
  full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
  save(full.res,file = paste0("~/output/Simulation_NG1000N200_homoSignal_V12345_mu",mu0,"_r20V1N270.RData"))
  
  res.tb.ls = lapply(full.res, "[[","res.tb")
  res.tb = data.frame(do.call(rbind,res.tb.ls))
  library(dplyr)
  df01 = res.tb %>% 
    group_by(orig_V,data.idx) %>%
    slice_max(V) %>%
    slice(which.max(avgR2_selected_soft_sepV))
  
  r2_cut = c(seq(0.01,0.1,0.01),seq(0.12,0.5,0.02))
  
  r2_cut_thrshld = data.frame(t(sapply(1:nrow(df01), function(i){
    res = full.res[[df01$l[i]]][["res"]]
    std.data = data.list[[df01$data.idx[i]]]
    pred_Glb = apply(res$postGV,1,function(x) {
      if(all(x == 0)){
        0
      }else{
        which.max(x)
      }
    })
    gene.lb.ls = lapply(1:res$V, function(v){which(pred_Glb == v)})
    if(res$V==1){
      mu_GK_0 = apply(std.data,1,mean)
      sigma_GK_0 = apply(std.data,1,sd)
      NG = nrow(std.data)
      N = ncol(std.data)
      
      pred_GN0 = matrix(rep(mu_GK_0,N),nrow = NG,ncol = N)
      SSE_G0 = apply((std.data-pred_GN0)^2, 1, sum)
      mu_GK_v = res$mu_GK_V[[1]]
      postvK = res$w_NK_V[[1]]
      pred_GNv = mu_GK_v %*% t(postvK)
      SSE_GV = apply((std.data-pred_GNv)^2, 1, sum)
      r2 = as.matrix(1-SSE_GV/SSE_G0)[pred_Glb == 1,1]
      
      R2_V_ls = list(data.frame(sort_r2 = sort(r2), genes = names(sort(r2)), v = 1) )
    }else{
      R2_V_ls = lapply(1:res$V, function(v) {
        r2 = res$R2_V[pred_Glb == v,v]
        df = data.frame(sort_r2 = sort(r2), genes = names(sort(r2)), v = v) 
      })
    }
    
    
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
    ratio.df$mu0 = as.numeric(df01$mu0[i])
    ratio.df$orig_V = as.numeric(df01$orig_V[i])
    ratio.df$data.idx = as.numeric(df01$data.idx[i])
    
    ratio.df_thrshld2 = ratio.df[which.min(abs(ratio.df$ratio-2)),]
    return(ratio.df_thrshld2)
  })))
  r2_cut_thrshld$mu0 = as.numeric(r2_cut_thrshld$mu0)
  r2_cut_thrshld$orig_V = as.numeric(r2_cut_thrshld$orig_V)
  r2_cut_thrshld$data.idx = as.numeric(r2_cut_thrshld$data.idx)
  r2cut.tab = r2_cut_thrshld
  save(r2cut.tab,
       file = paste0("~/output/Simulation_NG1000N200_homoSignal_V12345_mu",mu0,"_r2cutV1N270.RData"))
  
}


rm(list = ls())
options(stringsAsFactors = F)
library(Rcpp)
library(S4)
library(tightClust)
library(caret)
library(cluster)
library(MASS)
library(parallel)
sourceCpp("~/code/EM_MFClust_func_C.cpp")
source('~/code/EM_MFClust_func.R')
source('~/code/S4_Kmedoids.R')
source('~/code/run_func_selectR2.R')

V.list = c(1:2)
initial.list = expand.grid(seq(20,60,20),1:5)
colnames(initial.list) = c("kmin","seednum")
mu0 = c(0.8)
for (mu0 in c(0.8,1,1.2)) {
  load(paste0("~/data/Simulation_NG1000N200_homoSignal_datalist_",mu0,".RData"))
  load(paste0("~/output/Simulation_NG1000N200_homoSignal_V12345_mu",mu0,"_r2cutV1N270.RData"))
  
  set.seed(12345)
  scenario.idx = expand.grid(1:nrow(initial.list), 1:length(data.list),1:length(V.list))
  full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
  save(full.res,file = paste0("~/output/Simulation_NG1000N200_homoSignal_V12345_mu",mu0,"_selectR2V1N270.RData"))
  
}

