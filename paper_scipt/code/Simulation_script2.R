### Simulate data ##
rm(list=ls())
NG = 1000
N = 200
mu=1
p.list = seq(0,0.5,0.1)
for (i in 1:length(p.list)) {
  data.list = clust.lb0 = list()
  
  for(l in 1:20){
    set.seed(l)
    Y = matrix(rnorm(N*NG,0,1),nrow = NG,ncol = N)
    
    cl1 = sample(1:3,N,replace = T)
    cl2 = sample(1:3,N,replace = T)
    cl3 = sample(1:3,N,replace = T)
    
    m1 = m2 = m3 = mu
    
    #View1: 3clusters; 
    Y[1:30,cl1==2] = rnorm(30*sum(cl1==2),m1,2); Y[1:30,cl1==3] = rnorm(30*sum(cl1==3),-m1,2)
    Y[31:60,cl1==1] = rnorm(30*sum(cl1==1),m1,2); Y[31:60,cl1==2] = rnorm(30*sum(cl1==2),-m1,2)
    Y[61:90,cl1==3] = rnorm(30*sum(cl1==3),m1,2); Y[61:90,cl1==1] = rnorm(30*sum(cl1==1),-m1,2)
    
    #View2: 3clusters; 
    Y[91:120,cl2==2] = rnorm(30*sum(cl2==2),m2,2); Y[91:120,cl2==3] = rnorm(30*sum(cl2==3),-m2,2)
    Y[121:150,cl2==1] = rnorm(30*sum(cl2==1),m2,2); Y[121:150,cl2==2] = rnorm(30*sum(cl2==2),-m2,2)
    Y[151:180,cl2==3] = rnorm(30*sum(cl2==3),m2,2); Y[151:180,cl2==1] = rnorm(30*sum(cl2==1),-m2,2)
    
    #View3: 3clusters; 
    Y[181:210,cl3==2] = rnorm(30*sum(cl3==2),m3,2); Y[181:210,cl3==3] = rnorm(30*sum(cl3==3),-m3,2)
    Y[211:240,cl3==1] = rnorm(30*sum(cl3==1),m3,2); Y[211:240,cl3==2] = rnorm(30*sum(cl3==2),-m3,2)
    Y[241:270,cl3==3] = rnorm(30*sum(cl3==3),m3,2); Y[241:270,cl3==1] = rnorm(30*sum(cl3==1),-m3,2)
    
    #Fuzz genes
    Y[271:300,] = p.list[i]*Y[91:120,]+(1-p.list[i])*Y[181:210,]
    
    row.names(Y) = paste0("G",1:NG)
    colnames(Y) = paste0("S",1:N)
    
    data.list[[l]] = t(scale(t(Y)))
    
    
    clust.lb0[[l]] = list(v1 = cl1,
                          v2 = cl2,
                          v3 = cl3)
    gene.lb01 = list(v1 = 1:90, v2=c(91:180,271:300), v3 = c(181:270))
    gene.lb02 = list(v1 = 1:90, v2=c(91:180), v3 = c(181:270,271:300))
    
    gene.lb0_multi1 = c(rep(1,90),rep(2,90),rep(3,90),rep(2,30),rep(0,NG-300))
    gene.lb0_multi2 = c(rep(1,90),rep(2,90),rep(3,90),rep(2,30),rep(0,NG-300))
  }
  save(data.list,clust.lb0,gene.lb01,gene.lb02,
       file = paste0("~/data/Simulation_NG1000N200_homoSignal_datalist_p",p.list[i],"_fuzzmu1.RData"))
  
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
source('~/code/run_func_selectR2_fuzzy.R')

V.list = c(1:5)
initial.list = expand.grid(seq(20,60,20),1:5)
colnames(initial.list) = c("kmin","seednum")
for(p in seq(0,0.5,0.1)){
  mu0 = p #mu0 fixed, this is to replace the mu0 label in run script
  load(paste0("~/data/Simulation_NG1000N200_homoSignal_datalist_p",p,"_fuzzmu1.RData"))
  load(paste0("~/output/Simulation_NG1000N200_homoSignal_V12345_mu1_r2cut_NGv5.RData"))
  
  set.seed(12345)
  scenario.idx = expand.grid(1:nrow(initial.list), 1:length(data.list),1:length(V.list))
  full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
  save(full.res,file = paste0("~/output/Simulation_NG1000N200_homoSignal_V12345_p",mu0,"_selectR2_fuzzmu1.RData"))
  
}

