library(mclust)
EM_multiView = function(V, K, pV, pVK, mu_GK_V, sigma_GV, std.data, 
                        GCthres = 0.5, GCthres_num = 5, R2_cutoff = 0.2,
                        quite=F, updateK = T, updateK_thin = 1, maxr = 200, S4 = T){
  
  mu_GK_0 = apply(std.data,1,mean)
  sigma_GK_0 = apply(std.data,1,sd)
  
  NG = nrow(std.data)
  N = ncol(std.data)
  #:: parameters to update; 
  #:: V, K are decided internally by over-specify&thresholding and BIC respectively
  criterion = c()
  criterion.pV = c()
  criterion.pVK = c()
  criterion.mu_GK_V = c()
  criterion.sigma_GV = c()
  Kpath = list(K)
  delta_wgv = 10
  r = 1
  dis = 500
  # wGV.ls = list()
  # muGV.ls = list()
  # sigmaGV.ls = list()
  w_GV_old = 0
  while ((r<maxr & dis > 1e-3 & V>1)) {
    set.seed(r)
    if(quite == F){
      print(paste0("Run EM ",r," times; Distance ",dis))
    }
    #==E-STEP==#
    
    ### 1. w_GV ###
    
    ## gene specific cluster assignment probability
    
    logLik_VG = sapply(1:V, function(v){
      K_v = K[v]
      mu_GK = mu_GK_V[[v]]
      sigma_Gv = sigma_GV[,v]
      pK = pVK[[v]]
      pv = pV[v]
      logLik_vG2 = logLik_vG_hard2(N, NG, K_v, pv, std.data, mu_GK, sigma_Gv)
      
      # logLik_vG = sapply(1:nrow(mu_GK), function(g){
      #   mu_gK = mu_GK[g,]
      #   sigma_g = sigma_Gv[g]
      #   logLik_vgN = sum(sapply(1:N, function(n){
      #     xx = dnorm(std.data[g,n],mu_gK,sigma_g)
      #     log(xx)[which.max(xx)]
      #   }))
      #   #
      #   # t(sapply(1:N, function(n){
      #   #   xx = pK*dnorm(std.data[g,n],mu_gK,sigma_g)
      #   #   t_gn = xx/sum(xx)
      #   # }))
      #   #
      #   # t(sapply(1:N, function(n){
      #   #   xx = pK*dnorm(std.data[g,n],mu_gK,sigma_g)
      #   #   log(xx)
      #   # }))
      #   
      #   log(pV[v]) + logLik_vgN
      # })
      
      return(logLik_vG2)
    })
    
    
    w_GV = t(apply(logLik_VG,1,function(xx) {
      logScaleNorm = xx - max(xx)
      rmlogScaleNorm = exp(logScaleNorm)
      rmlogScaleNorm/sum(rmlogScaleNorm)
    }))
    
    # ### 2. w_NK_V ###
    #:: v from 1 to V i.e. V
    logScaleConditionalLikl_ls = ConditionalLL_ls_func(N, NG, V, K, std.data,
                                                       mu_GK_V, sigma_GV)
    w_NK_V = lapply(1:V, function(v){
      w_Gv = w_GV[,v]
      logScaleConditionalLikl_v = logScaleConditionalLikl_ls[[v]]
      wlogScaleConditionalLikl_v = sapply(logScaleConditionalLikl_v,function(xx){
        apply(xx,2,function(yy) sum(yy*w_Gv))
      })
      w_NK_v = t(sapply(1:nrow(wlogScaleConditionalLikl_v), function(i) {
        logScale = log(pVK[[v]])+wlogScaleConditionalLikl_v[i,]
        if(all(logScale == (-Inf))){
          w_NK_v_i = rep(0,K[v])
          w_NK_v_i[which.max(logScale)] = 1
        }else{
          logScaleNorm = logScale - max(logScale)
          rmlogScaleNorm = exp(logScaleNorm)
          w_NK_v_i = rmlogScaleNorm/sum(rmlogScaleNorm)
        }
        return(w_NK_v_i)
      }))
      
      return(w_NK_v)
    })
    
    pred_GN0 = matrix(rep(mu_GK_0,N),nrow = NG,ncol = N)
    SSE_G0 = apply((std.data-pred_GN0)^2, 1, sum)
    
    SSE_GV = sapply(1:V,function(v){
      mu_GK_v = mu_GK_V[[v]]
      postvK = w_NK_V[[v]]
      pred_GNv = mu_GK_v %*% t(postvK)
      SSE_G = apply((std.data-pred_GNv)^2, 1, sum)
      return(SSE_G)
    })
    # minSSE = apply(SSE_GV,1,min)
    # SSE0_cut = 1-qchisq(1-typeI_v0/NG,N-1) #after adjustment, needs really strong effect to reject
    
    R2_V = 1-apply(SSE_GV,2,function(x) x/SSE_G0)
    maxR2 = apply(R2_V,1,max)
    nonclust_idx = which(maxR2 < R2_cutoff)
    w_GV[nonclust_idx,] = matrix(c(rep(0,length(nonclust_idx)*V)), nrow = length(nonclust_idx))
    
    clust_idx = which(maxR2 >= R2_cutoff)
    if(length(clust_idx) < 4){
      return(list(V=V, K=K, pV=pV, pVK=pVK, mu_GK_V=mu_GK_V,  sigma_GV =sigma_GV,
                  Kpath=Kpath,
                  criterion.pV = criterion.pV, criterion.pVK = criterion.pVK,postGV = w_GV, w_NK_V = w_NK_V,
                  R2_cutoff = R2_cutoff,
                  avgR2_selected_soft = NA, avgR2_selected_hardView = NA,avgR2_selected_hardViewClust=NA,
                  minR2_selected_soft = NA, minR2_selected_hardView = NA, minR2_selected_hardViewClust = NA,
                  avgR2_selected_soft_sepV = NA, avgR2_selected_hardView_sepV = NA,avgR2_selected_hardViewClust_sepV=NA,
                  minR2_selected_soft_sepV = NA, minR2_selected_hardView_sepV = NA, minR2_selected_hardViewClust_sepV = NA,
                  criterion.mu_GK_V = criterion.mu_GK_V, criterion.sigma_GV = criterion.sigma_GV,
                  criterion = criterion, err = 1))
    }
    predGV_clust = apply(w_GV[clust_idx,], 1, which.max)
    tb = table(predGV_clust)
    null.v = union(as.numeric(names(tb)[tb<4]),which(apply(w_GV,2,function(x) sum(x > 1e-10) == 0)))
    #null.v = which(apply(w_GV,2,function(x) sum(x != 0)<2))
    if(length(null.v) != 0){
      V = V-length(null.v)
      if(V == 0){
        print("Clusterable views are reduced to 0")
        return(list(V=V, K=K, pV=pV, pVK=pVK, mu_GK_V=mu_GK_V,  sigma_GV =sigma_GV,
                    Kpath=Kpath,
                    criterion.pV = criterion.pV, criterion.pVK = criterion.pVK,postGV = w_GV, w_NK_V = w_NK_V,
                    R2_cutoff = R2_cutoff,
                    avgR2_selected_soft = NA, avgR2_selected_hardView = NA,avgR2_selected_hardViewClust=NA,
                    minR2_selected_soft = NA, minR2_selected_hardView = NA, minR2_selected_hardViewClust = NA,
                    avgR2_selected_soft_sepV = NA, avgR2_selected_hardView_sepV = NA,avgR2_selected_hardViewClust_sepV=NA,
                    minR2_selected_soft_sepV = NA, minR2_selected_hardView_sepV = NA, minR2_selected_hardViewClust_sepV = NA,
                    criterion.mu_GK_V = criterion.mu_GK_V, criterion.sigma_GV = criterion.sigma_GV,
                    criterion = criterion, err = 1))
      } 
      K = K[-null.v]
      w_GV = w_GV[,-(null.v),drop = F]
      idx = apply(w_GV,1,function(x){!all(x == 0)})
      if(ncol(w_GV) == 1){
        w_GV[idx,] = apply(w_GV[idx,,drop = F], 1, function(xx) xx/sum(xx))
      }else{
        w_GV[idx,] = t(apply(w_GV[idx,], 1, function(xx) xx/sum(xx)))
      }
      w_NK_V = w_NK_V[-null.v]
      mu_GK_V = mu_GK_V[-null.v]
      sigma_GV = sigma_GV[,-null.v,drop = F]
      pV = pV[-(null.v)]
      pV = pV/(sum(pV))
      pVK = pVK[-null.v]
      V = length(mu_GK_V)
      print(paste0("Adjust V to ", V,"!!!!!!!!!!!!!!!!!!!"))
    }
    
    if(updateK == T & delta_wgv > 0.1 & r %% updateK_thin == 0){
      for(v in 1:V){
        if(sum(w_GV[,v] > 1/V) < 5){#GCthres_num = 5
          subGmat = std.data[which(w_GV[,v] > 0),,drop=F]
        }else{
          subGmat = std.data[which(w_GV[,v] > 1/V),,drop=F]
        }
        if(K[v] == 2){
          Kcdt = c(K[v], K[v]+1)
        }else{
          Kcdt = c(K[v]-1, K[v], K[v]+1)
        }
        invisible(capture.output({
          res.S4 = K.Clust(t(subGmat),method="S4",Kmin=min(Kcdt),Kmax=max(Kcdt),trim.S4=0.05,cutoff=0,n.resample=50)
          aK_new = res.S4[[1]]
          res.Kmeans = kmeans(t(subGmat),centers=aK_new,nstart = 100)
          clust.lb = res.Kmeans$cluster
        }))
        
        if(all(table(res.Kmeans$cluster) != 1) & aK_new != K[v]){
          print(paste0("Update K in view ",v))
          prop = table(res.Kmeans$cluster)/sum(res.Kmeans$cluster)
          
          clust.fit.ls = lapply(1:aK_new, function(k){
            clust.expr = std.data[,clust.lb == k,drop = F]
            clust.mean = apply(clust.expr,1, function(m) {
              fitdistr(m, "normal")$estimate["mean"]
            })
            clust.sigma = apply(clust.expr,1, function(m) {
              fitdistr(m, "normal")$estimate["sd"]
            })
            return(list(clust.mean = clust.mean, clust.sigma = clust.sigma))
          })
          mean.mat = sapply(1:length(clust.fit.ls), function(x) clust.fit.ls[[x]]$clust.mean)
          sigma.mat = sapply(1:length(clust.fit.ls), function(x) clust.fit.ls[[x]]$clust.sigma)
          
          sigma.vec = apply(sigma.mat,1,function(x) sqrt(mean(x^2)))
          
          mu_GK_V[[v]] = mean.mat
          sigma_GV[,v] = sigma.vec
          pVK[[v]] = prop
          K[v] = aK_new
        }
        
      }
      Kpath = c(Kpath,list(K))
    }
    #reupdate cluster assignment
    logScaleConditionalLikl_ls = ConditionalLL_ls_func(N, NG, V, K, std.data,
                                                       mu_GK_V, sigma_GV)
    w_NK_V = lapply(1:V, function(v){
      w_Gv = w_GV[,v]
      logScaleConditionalLikl_v = logScaleConditionalLikl_ls[[v]]
      wlogScaleConditionalLikl_v = sapply(logScaleConditionalLikl_v,function(xx){
        apply(xx,2,function(yy) sum(yy*w_Gv))
      })
      w_NK_v = t(sapply(1:nrow(wlogScaleConditionalLikl_v), function(i) {
        logScale = log(pVK[[v]])+wlogScaleConditionalLikl_v[i,]
        if(all(logScale == (-Inf))){
          w_NK_v_i = rep(0,K[v])
          w_NK_v_i[which.max(logScale)] = 1
        }else{
          logScaleNorm = logScale - max(logScale)
          rmlogScaleNorm = exp(logScaleNorm)
          w_NK_v_i = rmlogScaleNorm/sum(rmlogScaleNorm)
        }
        return(w_NK_v_i)
      }))
      
      return(w_NK_v)
    })
    
    #==M-STEP==#
    
    ### 1. pV ###
    pV_new = colSums(w_GV)/sum(colSums(w_GV))
    
    ### 2. pVK ###
    pVK_new = lapply(1:V, function(v){
      colSums(w_NK_V[[v]])/N
    })
    
    ### 3. mu_g ###
    w_NKV = do.call(cbind,w_NK_V)
    mu_GKV = do.call(cbind,mu_GK_V)
    mu_GKV_new = t(t(std.data %*% w_NKV)/colSums(w_NKV))
    mu_GK_V_new = mat_to_list(m = mu_GKV_new, V=V, K=K)
    
    sigma_GV_new = sigma_GKV_func(mu_GK_V_new,std.data,w_NK_V,N,NG,V)
    
    
    #==check diff==#
    dis.pV = sum((as.numeric(pV_new)-as.numeric(pV))^2,na.rm = T)
    dis.pVK = sum((unlist(pVK_new)-unlist(pVK))^2,na.rm = T)
    dis.mu_GK_V = sum(sapply(1:V, function(v){
      Gindex = which(w_GV[,v] != 0)
      if(ncol(mu_GK_V_new[[v]]) == ncol(mu_GK_V_new[[v]]) & length(Gindex) != 0){
        adis = sum((mu_GK_V_new[[v]][Gindex,]-mu_GK_V[[v]][Gindex,])^2,na.rm =T)
      }else{
        adis = 500
      }
    }))
    dis.sigma_GV = sum((as.numeric(sigma_GV_new)-as.numeric(sigma_GV))^2,na.rm = T)
    dis = sqrt(dis.pV+dis.pVK+dis.mu_GK_V+dis.sigma_GV)
    
    criterion.pV[r] = dis.pV
    criterion.pVK[r] = dis.pVK
    criterion.mu_GK_V[r] = dis.mu_GK_V
    criterion.sigma_GV[r] = dis.mu_GK_V
    criterion[r] = dis
    
    #==record Kpath==#
    # Kpath = rbind(Kpath, K_new)
    delta_wgv = sum((as.numeric(w_GV_old)-as.numeric(w_GV))^2,na.rm = T)
    if(quite == F){
      print(paste0("wGV change: ",delta_wgv))
    }
    
    #==revalue==#
    # K = K_new
    pV = pV_new
    pVK = pVK_new
    mu_GK_V = mu_GK_V_new
    sigma_GV = sigma_GV_new
    w_GV_old = w_GV
    r = r+1
  }
  
  logLik_VG = sapply(1:V, function(v){
    K_v = K[v]
    mu_GK = mu_GK_V[[v]]
    sigma_Gv = sigma_GV[,v]
    pK = pVK[[v]]
    logLik_vG = sapply(1:nrow(mu_GK), function(g){
      mu_gK = mu_GK[g,]
      sigma_g = sigma_Gv[g]
      logLik_vgN = sum(sapply(1:N, function(n){
        xx = pK*dnorm(std.data[g,n],mu_gK,sigma_g)
        t_gn = xx/sum(xx)
        sum(t_gn*(log(xx)))
      }))
      log(pV[v]) + logLik_vgN
    })
    return(logLik_vG)
  })
  w_GV = t(apply(logLik_VG,1,function(xx) {
    logScaleNorm = xx - max(xx)
    rmlogScaleNorm = exp(logScaleNorm)
    rmlogScaleNorm/sum(rmlogScaleNorm)
  }))
  
  logScaleConditionalLikl_ls = ConditionalLL_ls_func(N, NG, V, K, std.data,
                                                     mu_GK_V, sigma_GV)
  w_NK_V = lapply(1:V, function(v){
    w_Gv = w_GV[,v]
    logScaleConditionalLikl_v = logScaleConditionalLikl_ls[[v]]
    wlogScaleConditionalLikl_v = sapply(logScaleConditionalLikl_v,function(xx){
      apply(xx,2,function(yy) sum(yy*w_Gv))
    })
    w_NK_v = t(sapply(1:nrow(wlogScaleConditionalLikl_v), function(i) {
      logScale = log(pVK[[v]])+wlogScaleConditionalLikl_v[i,]
      if(all(logScale == (-Inf))){
        w_NK_v_i = rep(0,K[v])
        w_NK_v_i[which.max(logScale)] = 1
      }else{
        logScaleNorm = logScale - max(logScale)
        rmlogScaleNorm = exp(logScaleNorm)
        w_NK_v_i = rmlogScaleNorm/sum(rmlogScaleNorm)
      }
      return(w_NK_v_i)
    }))
    
    return(w_NK_v)
  })
  
  pred_GN0 = matrix(rep(mu_GK_0,N),nrow = NG,ncol = N)
  SSE_G0 = apply((std.data-pred_GN0)^2, 1, sum)
  
  SSE_GV = sapply(1:V,function(v){
    mu_GK_v = mu_GK_V[[v]]
    postvK = w_NK_V[[v]]
    pred_GNv = mu_GK_v %*% t(postvK)
    SSE_G = apply((std.data-pred_GNv)^2, 1, sum)
    return(SSE_G)
  })
  
  SSE_GV_hard = sapply(1:V,function(v){
    mu_GK_v = mu_GK_V[[v]]
    postvK = w_NK_V[[v]]
    postvK_hard = apply(postvK, 1, function(xx) {
      out = rep(0,length(xx))
      out[which.max(xx)] = 1
      return(out)
    })
    pred_GNv = mu_GK_v %*% postvK_hard
    SSE_G = apply((std.data-pred_GNv)^2, 1, sum)
    return(SSE_G)
  })
  
  if(V != 1){
    R2_V = 1-apply(SSE_GV,2,function(x) x/SSE_G0)
    maxR2 = apply(R2_V,1,max)
    nonclust_idx = which(maxR2 < R2_cutoff)
    w_GV[nonclust_idx,] = 0
    clust_idx = which(maxR2 >= R2_cutoff)
    predGV_clust = apply(w_GV[clust_idx,], 1, which.max)
    
    R2_V_hard = 1-apply(SSE_GV_hard,2,function(x) x/SSE_G0)
    
  }else{
    R2_V = as.matrix(1-SSE_GV/SSE_G0)
    maxR2 = R2_V[,1]
    nonclust_idx = which(maxR2 < R2_cutoff)
    w_GV[nonclust_idx,] = 0
    clust_idx = which(maxR2 >= R2_cutoff)
    predGV_clust = rep(1,length(clust_idx))
    
    R2_V_hard = as.matrix(1-SSE_GV_hard/SSE_G0)
    
  }
  
  tb = table(predGV_clust)
  null.v = union(as.numeric(names(tb)[tb<4]),which(apply(w_GV,2,function(x) sum(x > 1e-10) == 0)))
  if(length(null.v) != 0){
    V = V-length(null.v)
    if(V == 0){
      print("Clusterable views are reduced to 0, consider a smaller R2 cutoff")
      return(list(V=V, K=K, pV=pV, pVK=pVK, mu_GK_V=mu_GK_V,  sigma_GV =sigma_GV,
                  Kpath=Kpath, 
                  criterion.pV = criterion.pV, criterion.pVK = criterion.pVK, postGV = w_GV, w_NK_V = w_NK_V,
                  R2_cutoff = R2_cutoff,
                  avgR2_selected_soft = NA, avgR2_selected_hardView = NA,avgR2_selected_hardViewClust=NA,
                  minR2_selected_soft = NA, minR2_selected_hardView = NA, minR2_selected_hardViewClust = NA,
                  avgR2_selected_soft_sepV = NA, avgR2_selected_hardView_sepV = NA,avgR2_selected_hardViewClust_sepV=NA,
                  minR2_selected_soft_sepV = NA, minR2_selected_hardView_sepV = NA, minR2_selected_hardViewClust_sepV = NA,
                  criterion.mu_GK_V = criterion.mu_GK_V, criterion.sigma_GV = criterion.sigma_GV,
                  criterion = criterion,err = 1))
    } 
    K = K[-null.v]
    w_GV = w_GV[,-(null.v),drop = F]
    idx = apply(w_GV,1,function(x){!all(x == 0)})
    if(ncol(w_GV) == 1){
      w_GV[idx,] = apply(w_GV[idx,,drop = F], 1, function(xx) xx/sum(xx))
    }else{
      w_GV[idx,] = t(apply(w_GV[idx,], 1, function(xx) xx/sum(xx)))
    }    
    w_NK_V = w_NK_V[-null.v]
    mu_GK_V = mu_GK_V[-null.v]
    sigma_GV = sigma_GV[,-null.v,drop = F]
    pV = pV[-(null.v)]
    pV = pV/(sum(pV))
    pVK = pVK[-null.v]
    R2_V = R2_V[,-null.v,drop = F]
    R2_V_hard = R2_V_hard[,-null.v,drop = F]
    
    V = length(mu_GK_V)
    print(paste0("Adjust V to ", V,"!!!!!!!!!!!!!!!!!!!"))
  }
  if(V != 1){
    w_GV_hard = t(apply(w_GV, 1, function(xx) {
      out = rep(0,length(xx))
      out[which.max(xx)] = 1
      return(out)
    }))
  }else{
    w_GV_hard = w_GV
  }
  
  avgR2_selected_soft = mean(rowSums(w_GV * R2_V)[clust_idx])
  avgR2_selected_soft_sepV = mean(sapply(1:V, function(v){
    xx = which(apply(w_GV, 1, which.max) == v) 
    mean(rowSums(w_GV * R2_V)[intersect(xx,clust_idx)])
  }))
  
  avgR2_selected_hardView = mean(rowSums(w_GV_hard * R2_V)[clust_idx])
  avgR2_selected_hardView_sepV = mean(sapply(1:V, function(v){
    xx = which(apply(w_GV, 1, which.max) == v) 
    mean(rowSums(w_GV_hard * R2_V)[intersect(xx,clust_idx)])
  }))
  
  avgR2_selected_hardViewClust = mean(rowSums(w_GV_hard * R2_V_hard)[clust_idx])
  avgR2_selected_hardViewClust_sepV = mean(sapply(1:V, function(v){
    xx = which(apply(w_GV, 1, which.max) == v) 
    mean(rowSums(w_GV_hard * R2_V_hard)[intersect(xx,clust_idx)])
  }))
  
  minR2_selected_soft = min(rowSums(w_GV * R2_V)[clust_idx])
  minR2_selected_soft_sepV = mean(sapply(1:V, function(v){
    xx = which(apply(w_GV, 1, which.max) == v) 
    min(rowSums(w_GV * R2_V)[intersect(xx,clust_idx)])
  }))
  
  minR2_selected_hardView = min(rowSums(w_GV_hard * R2_V)[clust_idx])
  minR2_selected_hardView_sepV = mean(sapply(1:V, function(v){
    xx = which(apply(w_GV, 1, which.max) == v) 
    min(rowSums(w_GV_hard * R2_V)[intersect(xx,clust_idx)])
  }))
  
  minR2_selected_hardViewClust = min(rowSums(w_GV_hard * R2_V_hard)[clust_idx])
  minR2_selected_hardViewClust_sepV = mean(sapply(1:V, function(v){
    xx = which(apply(w_GV, 1, which.max) == v) 
    min(rowSums(w_GV_hard * R2_V_hard)[intersect(xx,clust_idx)])
  }))
  return(list(V=V, K=K, pV=pV, pVK=pVK, mu_GK_V=mu_GK_V,  sigma_GV =sigma_GV,
              Kpath=Kpath, 
              criterion.pV = criterion.pV, criterion.pVK = criterion.pVK,
              criterion.mu_GK_V = criterion.mu_GK_V, criterion.sigma_GV = criterion.sigma_GV,
              criterion = criterion, R2_V = R2_V, postGV = w_GV, w_NK_V = w_NK_V, R2_cutoff = R2_cutoff,
              avgR2_selected_soft = avgR2_selected_soft, avgR2_selected_hardView = avgR2_selected_hardView,avgR2_selected_hardViewClust=avgR2_selected_hardViewClust,
              minR2_selected_soft = minR2_selected_soft, minR2_selected_hardView = minR2_selected_hardView, minR2_selected_hardViewClust = minR2_selected_hardViewClust,
              avgR2_selected_soft_sepV = avgR2_selected_soft_sepV, avgR2_selected_hardView_sepV = avgR2_selected_hardView_sepV,avgR2_selected_hardViewClust_sepV=avgR2_selected_hardViewClust_sepV,
              minR2_selected_soft_sepV = minR2_selected_soft_sepV, minR2_selected_hardView_sepV = minR2_selected_hardView_sepV, minR2_selected_hardViewClust_sepV = minR2_selected_hardViewClust_sepV,
              err = 0))#,wGV.ls = wGV.ls, muGV.ls = muGV.ls, sigmaGV.ls = sigmaGV.ls))
}


