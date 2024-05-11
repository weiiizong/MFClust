#1. V = 1:7, eta(R2cut)=0 =================================================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(tightClust)
library(cluster)
library(MASS)

load("~/data/caudate_putamen_PreprocessedData.RData")
sourceCpp("~/code/EM_MFClust_func_C.cpp")
source('~/code/EM_MFClust_func.R')
source('~/code/S4_Kmedoids.R')
source('~/code/run_func_initKmedoids_DSpooleddata_V.R')

std.data = std.data2000
V.list = 1:7
R2.list = c(0)
G.K.list = c(20)
initial.list = expand.grid(c(50,100,500),1:10)
colnames(initial.list) = c("kmin","seednum")
scenario.idx = expand.grid(1:nrow(initial.list), 1:length(R2.list),
                           1:length(G.K.list),1:length(V.list))
library(parallel)
full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
save(full.res,file = "~/output/caudate_putamen_pooled_full.res_V234567_r20.RData")


rm(list = ls())
load("~/data/caudate_putamen_PreprocessedData.RData")
std.data = std.data2000
load("~/output/caudate_putamen_pooled_full.res_V234567_r20.RData")

res.tb.ls = lapply(full.res, "[[","res.tb")
res.tb = data.frame(do.call(rbind,res.tb.ls))
library(dplyr)
df01 = res.tb %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice_max(V) %>%
  slice_max(avgR2_selected_soft_sepV)

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
     file = "~/output/caudate_putamen_full.res_V234567_avg_r2_cut_r20_2folds.RData") 

#2. V = 1:7, selected R2cut =================================================================================================
rm(list=ls())
library(Rcpp)
library(S4)
library(tightClust)
library(cluster)
library(MASS)

load("~/data/caudate_putamen_PreprocessedData.RData")
sourceCpp("~/code/EM_MFClust_func_C.cpp")
source('~/code/EM_MFClust_func.R')
source('~/code/S4_Kmedoids.R')
source('~/code/run_func_initKmedoids_DSpooleddata_V_selectR2_cutoff.R')

load("~/output/caudate_putamen_full.res_V234567_avg_r2_cut_r20_2folds.RData")

R2.list = avg_r2_cut
std.data = std.data2000
V.list = 1:7
G.K.list = c(20)
initial.list = expand.grid(c(50,100,500),1:10)
colnames(initial.list) = c("kmin","seednum")
scenario.idx = expand.grid(1:nrow(initial.list), 1:length(G.K.list),1:length(V.list))
library(parallel)
full.res = mclapply(1:nrow(scenario.idx),run.func_V,mc.cores = 50)
save(full.res,file = "~/output/caudate_putamen_full.res_V234567_selectR2_cutoff_2folds.RData")

#3. Plot ========================================================================================================================
load("~data/caudate_putamen_PreprocessedData.RData")
load("~output/caudate_putamen_full.res_V234567_selectR2_cutoff_2folds.RData")
std.data = std.data2000
res.tb.ls = lapply(full.res, "[[","res.tb")
res.tb = data.frame(do.call(rbind,res.tb.ls))
library(dplyr)
df01 = res.tb %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice_max(V) %>%
  slice_max(avgR2_selected_soft_sepV)
df01 %>% dplyr::select(c("orig_V","V","initial.kmin","max_pairARI","avg_pairARI","avgR2_selected_soft_sepV"))
res = full.res[[df01$l[5]]][["res"]]

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
gene.ls.final = gene.ls
cl.ls = lapply(1:res$V, function(v){
  apply(res$w_NK_V[[v]],1,which.max)
})
cl.ls.final = cl.ls

## Heatmap
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
order_varls <- list(f1 = c("RIN"), f2 = c("Psychotic"),
                    f3 = c("Brain.Region"), f4 = c("pH","RIN"), f5=c("Age","Sex"),
                    f6= c("Brain.Region","pH","RIN"))
p.ls = list()
light_green <- rgb(144, 238, 144, maxColorValue=255)  # Light Green
dark_green <- rgb(0, 100, 0, maxColorValue=255)
for(v in 1:res$V){
  dat = std.data[gene.ls[[v]],]
  order_var <- order_varls[[v]]
  col_pheno = cbind(cl.ls[[v]], clinical[,order_var])
  colnames(col_pheno) <- c("Cluster", order_var)
  row.names(col_pheno) = colnames(dat)
  col_pheno <- as.data.frame(col_pheno)
  col_pheno$Cluster <- factor(cl.ls[[v]], levels = c(1,2,3,4), labels = c(1,2,3,4))
  if("Psychotic" %in% order_var){
    col_pheno$Psychotic = factor(col_pheno$Psychotic)
  }
  if("Sex" %in% order_var){
    col_pheno$Sex = factor(col_pheno$Sex)
  }
  if("Brain.Region" %in% order_var){
    col_pheno$Brain.Region = factor(col_pheno$Brain.Region)
  }
  annoCol = list("Cluster"=setNames(c("grey","#D82390", "#E9D2E0", "#5E22DF"), levels(col_pheno$Cluster)),
                 "Psychotic"=setNames(c("red", "yellow"), levels(col_pheno$Psychotic)),
                 "Sex"=setNames(c("#1A6DF3", "#C5D8F7"), levels(col_pheno$Sex)),
                 "RIN" = colorRampPalette(c(dark_green, light_green))(100),
                 "pH"= colorRampPalette(c("red", "orange", "yellow"))(100),
                 "Age" = colorRampPalette(c("blue", "cyan"))(100),
                 "Brain.Region" = setNames(c("#5E22DF", "#CEC6E1"), levels(col_pheno$Brain.Region)))
  
  #HC within each cluster to determine order
  cl.ordered = list()
  for(i in 1:length(unique(cl.ls[[v]]))){
    dati = dat[,cl.ls[[v]] == i]
    xx = hclust(dist(t(dati)))
    cl.ordered[[i]] =  colnames(dati)[xx$order]
  }
  
  dat = apply(dat, 1, standardize)
  
  png(paste0("~/output/caudate_heatmap_v",v,".png"),width = 900,height = 850,
      res=200)
  pheatmap(t(dat)[,unlist(cl.ordered)], cluster_rows = T,cluster_cols = F, show_rownames = F, show_colnames = F,border_color = NA,
           color=colorRampPalette(c("purple","blue","yellow","orange"))(n = 200),treeheight_row = 0, treeheight_col = 0, labels_col="Samples",
           annotation_col = col_pheno,annotation_colors = annoCol,main = paste0("Facet ",v," (NG=",length(gene.ls[[v]]),")"),
           annotation_legend = T, annotation_names_col = T,legend = F)
  dev.off()
}
dev.off()

## Transition plot
df01_1 = res.tb1[which.max(res.tb1$avgR2_selected_soft_sepV),]

res.tb.ls2 = lapply(full.res, "[[","res.tb")
res.tb2 = data.frame(do.call(rbind,res.tb.ls2))
df01_2 = res.tb2 %>% filter(max_pairARI<0.2) %>%
  group_by(orig_V) %>%
  slice_max(V) %>%
  slice_max(avgR2_selected_soft_sepV)

clinical <- as.data.frame(clinical)
fac_var <- c("Brain.Region","Sex","Psychotic")
num_var <- c("Age","pH","RIN")

library("RColorBrewer")
library(pheatmap)
p.ls = list()
p.ls.V = cl.ls.V = gene.ls.V = list()
for (l in 1:7) {
  if(l == 1){
    s.idx = df01_1$l
    res = full.res_V1[[s.idx]][["res"]]
    V = full.res_V1[[s.idx]][["res.tb"]]["orig_V"]
  }else{
    s.idx = df01_2$l[l-1]
    res = full.res[[s.idx]][["res"]]
    V = full.res[[s.idx]][["res.tb"]]["orig_V"]
  }
  
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
  gene.ls.V[[l]] =  gene.ls
  cl.ls = lapply(1:res$V, function(v){
    apply(res$w_NK_V[[v]],1,which.max)
  })
  
  assoc <- list()
  for (var in fac_var) {
    clin <- clinical[,var]
    valid.id <- which(!is.na(clin))
    clin <- clin[valid.id]
    pv.ls <- lapply(1:length(cl.ls), function(v){
      cl <- as.factor(cl.ls[[v]])
      cl <- cl[valid.id]
      if(l==1){
        pv <- fisher.test(clin, cl, simulate.p.value=TRUE, B=1e6)$p.value
      }else{
        pv <- fisher.test(clin, cl)$p.value
      }
      return(pv)
    })
    assoc[[var]] <- unlist(pv.ls)
  }
  
  for (var in num_var) {
    clin <- as.numeric(clinical[,var])
    pv.ls <- lapply(1:length(cl.ls), function(v){
      cl <- as.factor(cl.ls[[v]])
      pv <- kruskal.test(clin ~ cl)$p.value
      return(pv)
    })
    assoc[[var]] <- unlist(pv.ls)
  }
  variableP <- do.call(rbind, assoc)
  colnames(variableP) <- paste0("View_", 1:res$V)
  if(l == 6){
    metric.tb = variableP
  }
  
  metric.tblog10 <- -log10(variableP)
  metric.tblog10 <- as.data.frame(metric.tblog10)
  metric.tblog10$variables <- rownames(metric.tblog10)
  p.ls = list()
  for(v in 1:res$V){
    df <- metric.tblog10[,c(paste0("View_", v), "variables")]
    colnames(df)[1] <- "log10P"
    if(((l==5)&(v==3))|((l %in% c(6,7))&(v==5))){
      p.ls[[v]] = ggplot(df,aes(x = variables, y = log10P, fill = variables))+
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 2.86, color = "red")+
        theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust=1, margin = margin(t = -20, unit = "pt")))+
        theme(axis.title.y = element_blank())+
        labs(title = paste0("Facet ",v))+
        scale_x_discrete(name = "")+
        theme(plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "none",
              axis.text=element_text(size=18),axis.title=element_text(size=14))+
        labs(fill = expression(paste("Cluster label association (", -log[10](P),")")))
    }else{
      p.ls[[v]] = ggplot(df,aes(x = variables, y = log10P, fill = variables))+
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 2.86, color = "red")+
        theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust=1, margin = margin(t = -20, unit = "pt")))+
        theme(axis.title.y = element_blank())+
        ylim(c(0,21))+
        labs(title = paste0("Facet ",v))+
        scale_x_discrete(name = "")+
        theme(plot.title = element_text(hjust = 0.5))+
        theme(legend.position = "none",
              axis.text=element_text(size=18),axis.title=element_text(size=14))+
        labs(fill = expression(paste("Cluster label association (", -log[10](P),")")))
    }
  }
  p.ls.V[[l]] <- p.ls
  cl.ls.V[[l]] = cl.ls
}
p.ls.V

## Barplot
grid.arrange(arrangeGrob(grobs = p.ls.V[[6]], ncol=3))

## Boxplot
variables <- rownames(metric.tb)
gene.ls = gene.ls.final
logp.ls = lapply(1:6, function(v){
  sub.std.data = std.data[gene.ls[[v]],,drop=F]
  log.df = data.frame(t(apply(sub.std.data, 1, function(y){
    sapply(variables, function(var){
      -log10(summary(lm(y ~ clinical[,var,drop = T]))$coefficients[2,4])
    })
  })))
  log.df$v = v
  return(log.df)
})
logp.df = do.call(rbind,logp.ls)

p.ls = list()
library(ggpubr)
for(v in 1:6){
  metric.tb.v = metric.tb[,v]
  names(metric.tb.v) = rownames(metric.tb)
  p.df = logp.df[logp.df$v == v,]
  df = gather(p.df,variable,log10P,-v)
  df$sigVar = sapply(1:nrow(df), function(xx) -log10(metric.tb.v[df$variable[xx]]))
  df$sigVar[df$sigVar>6] = 6
  ylm = ifelse(v == 5, 320, 40)
  p.ls[[v]] = ggplot(df,aes(x = variable, y = log10P,fill = sigVar))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title = paste0("Facet ",v," (NG=",length(gene.ls[[v]]),")"))+
    scale_fill_gradient(
      limits = c(0,6),
      low = "#C7D1D3",
      high = "#18C3E6",
      space = "Lab",
      na.value = "#C7D1D3",
      guide = "colourbar",
      aesthetics = "fill"
    )+
    ylim(c(0,ylm))+
    scale_x_discrete(name = "")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position = "none",axis.text=element_text(size=11),axis.title=element_text(size=13),
          axis.title.y = element_blank())
}
grid.arrange(arrangeGrob(grobs= p.ls,ncol=3))
