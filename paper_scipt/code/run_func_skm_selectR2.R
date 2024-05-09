run_func = function(m){
  data.idx = m
  std.data = data.list[[data.idx]]
  aclust.lb0 = clust.lb0[[data.idx]]
  
  R2_cutoff = r2avgV$r_avg[r2avgV$data.idx == data.idx]
  
  out = tryCatch({  
    data = std.data
    k_vector = 2:6
    skm.ls = list()
    
    for (v in 1:V) { 
      set.seed(v)
      wbounds_list = list()
      for(l in 1:length(k_vector)){#for each K, using the algorithm to get 20 lambda.
        wbounds_list[[l]] = region.lambda(lam1=1.1,iteration=20,t(data),k_vector[l])
      }
      for(l in 1:length(k_vector)){
        temp<-KMeansSparseCluster(t(data),K=k_vector[l],wbounds=wbounds_list[[l]],nstart=100)
        num<-rep(0,length(temp))
        for(i in 1:length(num)){#get the corresponding number of features for each K and each lambda
          num[i]<-sum(temp[[i]]$ws>0)
        }
        if(sum(num==ncol(data))>0){#For each K, if a certain lambda selects all features, delete it. 
          #For a large simulation study, two closest lambda next to it can also be removed to be conservative.
          wbounds_list[[l]]<-wbounds_list[[l]][1:(min(which(num==nrow(data)))-3)]
        }
      }
      res.S4<-KL.S4(x=t(data),lambda_list = wbounds_list,k_vector = k_vector,trim =0.05,n.resample = 50,num.cores = 1)
      
      km.out =  KMeansSparseCluster(t(data),K=res.S4$optimal_k,nstart = 150,silent = T)
      R2_check = lapply(1:length(km.out), function(i){
        gene.idx = km.out[[i]]$ws > 0
        adata = data[gene.idx,]
        pred.mat = adata
        for(k0 in 1:res.S4$optimal_k){
          mu = apply(adata[,which(km.out[[i]]$Cs==k0)],1,mean)
          pred.mat[,km.out[[i]]$Cs==k0] = matrix(c(rep(mu,sum(km.out[[i]]$Cs==k0))),nrow = nrow(adata),ncol = sum(km.out[[i]]$Cs==k0))
        }
        #check R2
        SSE_G = apply((pred.mat-adata)^2, 1, sum)
        pred_GN0 = matrix(rep(apply(adata,1,mean),ncol(adata)),nrow = nrow(adata),ncol = ncol(adata))
        SSE_G0 = apply((adata-pred_GN0)^2, 1, sum)
        minR2 = min(1 - SSE_G/SSE_G0)
        
        return(list(minR2 = minR2, NG = sum(gene.idx), pred.mat = pred.mat))
      })
      
      minR2 = sapply(R2_check, function(x) x$minR2)
      idx = which(minR2 > R2_cutoff)
      if(length(idx) == 0){
        next
      }else{
        select_idx = max(which(minR2>R2_cutoff))
        skm.ls[[v]] = km.out[[select_idx]]
        
        # pred.dat = data
        # for(k0 in 1:K[v]){
        #   mu = apply(data[,which(km.out[[select_idx]]$Cs==k0)],1,mean)
        #   pred.dat[,km.out[[select_idx]]$Cs==k0] = matrix(c(rep(mu,sum(km.out[[select_idx]]$Cs==k0))),nrow = nrow(data),ncol = sum(km.out[[select_idx]]$Cs==k0))
        # }
        # data = data - pred.dat
        
        amean.mat = sapply(1:res.S4$optimal_k, function(k){
          mu = apply(data[,which(km.out[[select_idx]]$Cs==k)],1,mean)
        })
        
        data = sapply(1:ncol(data), function(l){
          x = data[,l]
          cl = km.out[[select_idx]]$Cs[l]
          mu = amean.mat[,cl]
          x_new = x - mu*((t(mu) %*% x/t(mu) %*% mu)[1,1])
        })
      }    
    }
    if(length(skm.ls) == 0){
      aout = list(res.tb = NA,skm.ls = skm.ls, m = m)
    }else{
      len = sum(sapply(skm.ls, is.null))
      if(len == length(skm.ls)){
        aout = list(res.tb = NA,skm.ls = skm.ls, m = m)
      }else{
        if(len == 0){
          skm.ls = skm.ls
        }else{
          skm.ls = skm.ls[!sapply(skm.ls, is.null)]
        }
        library(mclust)
        pred.gene.lb = lapply(skm.ls,function(xx) which(xx$ws > 0))
        pred.cl.lb = lapply(skm.ls, function(xx) xx$Cs)
        res.tb = t(sapply(1:length(skm.ls), function(v){
          cl = pred.cl.lb[[v]]
          ari_res = sapply(1:length(aclust.lb0), function(xx){
            adjustedRandIndex(aclust.lb0[[xx]],cl)
          })
          Cl_res = c(ARI = max(ari_res), which = which.max(ari_res))
          orig.idx = which.max(ari_res)
          
          binary_lb = factor(ifelse(1:1000 %in% pred.gene.lb[[v]], 1, 0), levels = c("1","0"))
          true.gene.lb = rep(0, nrow(std.data))
          true.gene.lb[gene.lb0[[orig.idx]]] = 1
          true.gene.lb = factor(true.gene.lb, levels = c("1","0"))
          mat = confusionMatrix(binary_lb,true.gene.lb,positive = "1")
          G_Sensitivity = mat$byClass["Sensitivity"]
          G_Specificity = mat$byClass["Specificity"]
          
          return(c(mu0 = mu0,R2cutoff = R2_cutoff, v = v, orig.idx = orig.idx, Cl_ARI = max(ari_res), NG = length(pred.gene.lb[[v]]),
                   G_Sensitivity, G_Specificity,m = m,data.idx = data.idx))
        }))
        aout = list(res.tb = res.tb,skm.ls = skm.ls,m = m)
        
        
      }
    }
    aout
  },
  error = function(e){
    return(list(e = e,m = m))
  })
  return(out)
}
