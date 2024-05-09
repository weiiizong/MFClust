#skm ========================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(parallel)
library(sparcl)
library(caret)
library(dplyr)
source('~/code/run_func_skm_selectR2.R')

for(mu0 in c(0.8,1,1.2)){
  V = 5
  load(paste0("~/data/Simulation_NG1000N200_homoSignal_datalist_",mu0,".RData"))
  load(paste0("~/output/Simulation_NG1000N200_homoSignal_V12345_mu",mu0,"_r2cut.RData"))
  r2cut.tab$r = as.numeric(r2cut.tab$r)
  r2avgV = r2cut.tab%>%
    group_by(data.idx) %>%
    summarise(r_avg=mean(r))
  
  set.seed(12345)
  skm.res.ls = mclapply(1:length(data.list),run_func,mc.cores = 50)
  save(skm.res.ls,file = paste0("~/output/Simulation_NG1000N200_homoSignal_V5_mu_",mu0,"_selectR2_skm.RData"))
  
}
#sgm ========================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(parallel)
library(sparcl)
library(caret)
library(Brobdingnag)
library(dplyr)
source('~/code/run_func_sgm_selectR2.R')

V = 5
for(mu0 in c(0.8,1,1.2)){
  load(paste0("~/data/Simulation_NG1000N200_homoSignal_datalist_",mu0,".RData"))
  load(paste0("~/output/Simulation_NG1000N200_homoSignal_V12345_mu",mu0,"_r2cut.RData"))
  r2cut.tab$r = as.numeric(r2cut.tab$r)
  r2avgV = r2cut.tab%>%
    group_by(data.idx) %>%
    summarise(r_avg=mean(r))
  
  set.seed(12345)
  sgm.res.ls = mclapply(1:length(data.list),run_func,mc.cores = 50)
  save(sgm.res.ls,file = paste0("~/ouput/Simulation_NG1000N200_homoSignal_V5_mu_",mu0,"_selectR2_sgm.RData"))
}