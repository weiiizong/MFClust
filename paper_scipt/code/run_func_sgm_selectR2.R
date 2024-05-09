run_func = function(m){
  data.idx = m
  
  std.data = data.list[[data.idx]]
  aclust.lb0 = clust.lb0[[data.idx]]
  
  R2_cutoff = r2avgV$r_avg[r2avgV$data.idx == data.idx]
  
  out = tryCatch({  
    data = std.data
    k_vector = 2:6
    sgm.ls = list()
    
    for (v in 1:V) {
      print(v)
      set.seed(v)
      res.ls = lapply(k_vector,function(k){
        #prepare initials
        km.perm = KMeansSparseCluster.permute(t(data),K=k,nperms=20,silent = T)
        km.out =  KMeansSparseCluster(t(data),K=k,nstart = 150,silent = T)
        R2_check = lapply(1:length(km.out), function(i){
          gene.idx = km.out[[i]]$ws > 0
          adata = data[gene.idx,]
          pred.mat = adata
          for(k0 in 1:k){
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
        if(length(idx) != 0){
          idx = max(which(minR2 >= R2_cutoff))
        }else{
          idx = which.min(R2_cutoff-minR2)
        }
        center.ls = lapply(1:k, function(k0){
          apply(data[,which(km.out[[idx]]$Cs==k0)],1,mean)
        })
        center_gauss = do.call(cbind, center.ls)
        
        tuning_param_gauss = seq(0,10,1)
        model_gauss = lapply(1:length(tuning_param_gauss),function(i){
          res = sgClust(data=data,c_center=center_gauss,
                        lambda=tuning_param_gauss[i],K=k)
          return(res)
        })
        
        R2_check2 = lapply(1:length(model_gauss), function(i){
          print(i)
          gene.idx = which(apply(model_gauss[[i]]$result$mu,1,function(x) all(x != 0)))
          if(length(gene.idx) <=1){
            return(list(quantR2 = 0))
          }else{
            adata = data[gene.idx,]
            pred.mat = model_gauss[[i]]$result$mu[gene.idx,] %*% t(model_gauss[[i]]$result$z)
            #check R2
            SSE_G = apply((pred.mat-adata)^2, 1, sum)
            pred_GN0 = matrix(rep(apply(adata,1,mean),ncol(adata)),nrow = nrow(adata),ncol = ncol(adata))
            SSE_G0 = apply((adata-pred_GN0)^2, 1, sum)
            quantR2 = quantile(1 - SSE_G/SSE_G0, 0.05)
            return(list(quantR2 = quantR2, NG = length(gene.idx), pred.mat = pred.mat))
          }
        })
        quantR2 = sapply(R2_check2, function(x) x$quantR2)
        if(sum(quantR2 >= R2_cutoff) >0){
          select_idx = min(which(quantR2 >= R2_cutoff))
        }else{
          select_idx = which.min(R2_cutoff-quantR2)
        }
        return(model_gauss[[select_idx]])
      })
      
      #select K
      idx = which.min(sapply(res.ls,function(xx) xx$BIC))
      K = k_vector[[idx]]
      sgm.ls[[v]] = res.ls[[idx]]
      
      amean.mat = res.ls[[idx]]$result$mu
      data = sapply(1:ncol(data), function(l){
        x = data[,l]
        mu = amean.mat[,apply(res.ls[[idx]]$result$z, 1, which.max)[l]]
        x_new = x - mu*((t(mu) %*% x/t(mu) %*% mu)[1,1])
      })
    }
    if(length(sgm.ls) == 0){
      aout = list(res.tb = NA,sgm.ls = sgm.ls, m = m)
    }else{
      len = sum(sapply(sgm.ls, is.null))
      if(len == length(sgm.ls)){
        aout = list(res.tb = NA,sgm.ls = sgm.ls, m = m)
      }else{
        if(len == 0){
          sgm.ls = sgm.ls
        }else{
          sgm.ls = sgm.ls[!sapply(sgm.ls, is.null)]
        }
        library(mclust)
        pred.gene.lb = lapply(sgm.ls,function(xx) which(apply(xx$result$mu,1,function(x) all(x != 0))))
        pred.cl.lb = lapply(sgm.ls, function(xx) apply(xx$result$z,1,which.max))
        res.tb = t(sapply(1:length(sgm.ls), function(v){
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
                   G_Sensitivity, G_Specificity,m = m,data.idx=data.idx))
        }))
        aout = list(res.tb = res.tb,sgm.ls = sgm.ls,m = m)
        
        
      }
    }
    aout
  },
  error = function(e){
    return(list(e = e,m = m))
  })
  return(out)
}
