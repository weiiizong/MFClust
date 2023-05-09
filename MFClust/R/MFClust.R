#' Model-based multi-facet clustering algorithm
#'
#' Provide multiple clustering solutions to the data provided
#'
#' @param V the number of facets/clustering solutions.
#' @param K the initial cluster numbers.A vector of length V indicating the number of
#' clusters in each facet.
#' @param pV the initial prior facet assignment probability vector of the length V.
#' @param pVK a list of length V where each element is the initial prior cluster assignment
#' probability vector in v.
#' @param mu_GK_V a list of length V where each element is a matrix. Each element represents
#' the initial cluster center of the column cluster for the row gene/feature.
#' @param sigma_GV a matrix storing the initial cluster standard deviation where each row
#' corresponds to a gene/feature and each column corresponds to a facet.
#' @param std.data a standardized data matrix. Rows are the genes/features and columns are
#' samples. Each row is standardized to mean zero and standard deviation one.
#' @param R2_cutoff genes/features whose fitted R2 is smaller than this value will be labeled
#' as non-clusterable genes and do not belong to any facet.
#' @param quite F if report each EM algorithm run.
#' @param updateK T if update the number of clusters during the procedure.
#' @param updateK_thin the number of EM runs before each K update. Default is 1.
#' @param maxr the maximum number of EM runs. Default is 200.
#'
#' @return A list object
#' \itemize{
#'  \item{\code{V}}{V}
#'  \item{\code{K}}{fitted K}
#'  \item{\code{pV}}{fitted pV}
#'  \item{\code{pVK}}{fitted pVK}
#'  \item{\code{mu_GK_V}}{fitted mu_GK_V}
#'  \item{\code{sigma_GV}}{fitted sigma_GV}
#'  \item{\code{Kpath}}{a record of K changes}
#'  \item{\code{postGV}}{posterior facet assignment for every gene/feature}
#'  \item{\code{postNK_V}}{posterior cluster assignment for every sample in each facet}
#'  \item{\code{R2_cutoff}}{genes/features whose fitted R2 is smaller than this value will be labeled as non-clusterable genes and do not belong to any facet.}
#'  \item{\code{avgLLik_selected_soft_sepV}}{the average likelihood of selected genes in each facet}
#'  \item{\code{avgR2_selected_soft_sepV}} {the average R2 of selected genes in each facet}
#'  \item{\code{criterion}} {a record of parameter changes throughout the EM runs}
#' }
#' @export MFClust
#'
#' @examples
#' \dontrun{
#' data(BA11_BA47_NG2000)
#' init.ls = MFClust_init(std.data=BA11_BA47_NG2000, V=5, G.K=20, initial.kmin=100,
#' R2_cutoff = 0.22, seed = 9)
#' V = length(init.ls)
#' pV=rep(1/V,V)
#' K=sapply(init.ls[1:V], function(xx) ncol(xx$mean.mat))
#' pVK=lapply(1:V, function(v) rep(1/K[v],K[v]))
#' mu_GK_V=lapply(init.ls[1:V], function(xx) {
#' mat = xx$mean.mat
#' row.names(mat) = row.names(std.data)
#' return(mat)
#' })
#' sigma_GV=sapply(init.ls[1:V], function(xx) {
#'   mat = apply(xx$sigma.mat,1,function(x){sqrt(mean(x^2))})
#'   return(mat)
#' })
#' row.names(sigma_GV) = row.names(std.data)
#' res = MFClust(V, K, pV, pVK, mu_GK_V, sigma_GV, std.data, R2_cutoff = 0.22)
#' }
MFClust = function(V, K, pV, pVK, mu_GK_V, sigma_GV, std.data,
                   R2_cutoff = 0.2, quite=F, updateK = T,
                   updateK_thin = 1, maxr = 200){
  if(length(K) != V|length(pV) != V|length(pVK) != V|length(mu_GK_V) != V|ncol(sigma_GV) != V ){
    stop("The initial parameter dimensions are not matched to V")
  }
  mu_GK_0 = apply(std.data,1,mean)
  sigma_GK_0 = apply(std.data,1,sd)

  NG = nrow(std.data)
  N = ncol(std.data)
  criterion = c()
  criterion.pV = c()
  criterion.pVK = c()
  criterion.mu_GK_V = c()
  criterion.sigma_GV = c()
  Kpath = list(K)
  delta_wgv = 10
  r = 1
  dis = 500
  w_GV_old = 0
  while ((r<maxr & dis > 1e-3 & V>1)) {
    set.seed(r)
    if(quite == F){
      print(paste0("Run EM ",r," times; Distance ",dis))
    }
    #==E-STEP==#

    ### 1. w_GV_star ###
    ##:: GV level always consider v from 0 to V i.e. c(0,V)
    logScaleMarginalLikl_mat_ls = lapply(1:V, function(v){
      K_v = K[v]
      mu_GK = mu_GK_V[[v]]
      sigma_Gv = sigma_GV[,v]
      pK = pVK[[v]]
      MarginalLikl_mat = MarginalLL_v_func(N, NG, K_v, std.data,
                                           mu_GK, sigma_Gv, pK)
      return(MarginalLikl_mat)
    })
    LL_GVmat = sapply(logScaleMarginalLikl_mat_ls, rowSums)
    if(V == 1){
      w_GV_star = matrix(1,nrow = NG, ncol = V)
      row.names(w_GV_star) = row.names(std.data)
    }else{
      baseLL = sapply(1:V, function(v) {
        K_v = K[v]
        pvK = pVK[[v]]
        mu_GK_v = mu_GK_V[[v]]
        sigma_Gv = sigma_GV[,v]

        if(K_v>1){
          idx_pair = combn(K_v,2)
          diffmu2G = apply(idx_pair,2,function(xx) {
            a = xx[1]
            b = xx[2]
            (pvK[a] * pvK[b])*(mu_GK_v[,a]-mu_GK_v[,b])^2
          })
          scale_p = sum(apply(idx_pair,2,function(xx) {
            a = xx[1]
            b = xx[2]
            (pvK[a] * pvK[b])
          }))
          rowSums(diffmu2G)/(scale_p*2*sigma_Gv^2) + log(sigma_Gv*sqrt(2*pi))
        }else{
          rep(0,length(sigma_Gv))
        }
      })
      LL_GVmat2 = t(apply(LL_GVmat,1,function(x) log(pV)+x)) #adjust for prior
      LL_GVmat3 = LL_GVmat2 + N*baseLL #adjust for relative likelihood
      w_GV = t(apply(LL_GVmat3,1, function(x) {
        exp(x-max(x))/sum(exp(x-max(x)))
      }))
      w_GV_star = w_GV

    }

    # ### 2. w_NK_V ###
    #:: v from 1 to V i.e. V
    logScaleConditionalLikl_ls = ConditionalLL_ls_func(N, NG, V, K, std.data,
                                                       pVK, mu_GK_V, sigma_GV)
    w_NK_V = lapply(1:V, function(v){
      w_Gv = w_GV_star[,v]
      logScaleConditionalLikl_v = logScaleConditionalLikl_ls[[v]]
      wlogScaleConditionalLikl_v = sapply(logScaleConditionalLikl_v,function(xx){
        apply(xx,2,function(yy) sum(yy*w_Gv))
      })
      w_NK_v = t(sapply(1:nrow(wlogScaleConditionalLikl_v), function(i) {
        logScale = wlogScaleConditionalLikl_v[i,]
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

    R2_V = 1-apply(SSE_GV,2,function(x) x/SSE_G0)
    maxR2 = apply(R2_V,1,max)
    nonclust_idx = which(maxR2 < R2_cutoff)
    w_GV_star[nonclust_idx,] = matrix(c(rep(0,length(nonclust_idx)*V)), nrow = length(nonclust_idx))

    clust_idx = which(maxR2 >= R2_cutoff)
    if(length(clust_idx) < 4){
      return(list(V=V, K=K, pV=pV, pVK=pVK, mu_GK_V=mu_GK_V,  sigma_GV =sigma_GV,
                  Kpath=Kpath,postGV = w_GV_star, postNK_V = w_NK_V,
                  R2_cutoff = R2_cutoff,
                  avgLLik_selected_soft_sepV=NA,avgR2_selected_soft_sepV = NA,
                  criterion = criterion))
    }
    predGV_clust = apply(w_GV_star[clust_idx,], 1, which.max)
    tb = table(predGV_clust)
    null.v = union(as.numeric(names(tb)[tb<4]),which(apply(w_GV_star,2,function(x) sum(x > 1e-10) == 0)))
    if(length(null.v) != 0){
      V = V-length(null.v)
      if(V == 0){
        print("Clusterable views are reduced to 0, consider a smaller R2 cutoff")
        return(list(V=V, K=K, pV=pV, pVK=pVK, mu_GK_V=mu_GK_V,  sigma_GV =sigma_GV,
                    Kpath=Kpath,postGV = w_GV_star, postNK_V = w_NK_V,
                    R2_cutoff = R2_cutoff,
                    avgLLik_selected_soft_sepV=NA,avgR2_selected_soft_sepV = NA,
                    criterion = criterion))      }
      K = K[-null.v]
      w_GV_star = w_GV_star[,-(null.v),drop = F]
      idx = apply(w_GV_star,1,function(x){!all(x == 0)})
      if(ncol(w_GV_star) == 1){
        w_GV_star[idx,] = apply(w_GV_star[idx,,drop = F], 1, function(xx) xx/sum(xx))
      }else{
        w_GV_star[idx,] = t(apply(w_GV_star[idx,], 1, function(xx) xx/sum(xx)))
      }
      w_NK_V = w_NK_V[-null.v]
      mu_GK_V = mu_GK_V[-null.v]
      sigma_GV = sigma_GV[,-null.v,drop = F]
      pV = pV[-(null.v)]
      pV = pV/(sum(pV))
      pVK = pVK[-null.v]
      V = length(mu_GK_V)
      print(paste0("Adjust V to ", V,"!"))
    }

    if(updateK == T & delta_wgv > 0.1 & r %% updateK_thin == 0){
      for(v in 1:V){
        if(sum(w_GV_star[,v] > 1/V) < 5){#For smaller views, use all genes
          subGmat = std.data[which(w_GV_star[,v] > 0),,drop=F]
        }else{
          subGmat = std.data[which(w_GV_star[,v] > 1/V),,drop=F]
        }
        if(K[v] == 2){
          Kcdt = c(K[v], K[v]+1)
        }else if(K[v] >= 5){
          Kcdt = c(K[v]-1, K[v])
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
                                                       pVK, mu_GK_V, sigma_GV)
    w_NK_V = lapply(1:V, function(v){
      w_Gv = w_GV_star[,v]
      logScaleConditionalLikl_v = logScaleConditionalLikl_ls[[v]]
      wlogScaleConditionalLikl_v = sapply(logScaleConditionalLikl_v,function(xx){
        apply(xx,2,function(yy) sum(yy*w_Gv))
      })
      w_NK_v = t(sapply(1:nrow(wlogScaleConditionalLikl_v), function(i) {
        logScale = wlogScaleConditionalLikl_v[i,]
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
    pV_new = colSums(w_GV_star)/sum(colSums(w_GV_star))

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
      Gindex = which(w_GV_star[,v] != 0)
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

    delta_wgv = sqrt(sum((as.numeric(w_GV_old)-as.numeric(w_GV_star))^2,na.rm = T))
    if(quite == F){
      print(paste0("wGV change: ",delta_wgv))
    }

    #==revalue==#
    pV = pV_new
    pVK = pVK_new
    mu_GK_V = mu_GK_V_new
    sigma_GV = sigma_GV_new
    w_GV_old = w_GV_star
    r = r+1
  }
  logScaleMarginalLikl_mat_ls = lapply(1:V, function(v){
    K_v = K[v]
    mu_GK = mu_GK_V[[v]]
    sigma_Gv = sigma_GV[,v]
    pK = pVK[[v]]
    MarginalLikl_mat = MarginalLL_v_func(N, NG, K_v, std.data,
                                         mu_GK, sigma_Gv, pK)
    return(MarginalLikl_mat)
  })
  LL_GVmat = as.matrix(sapply(logScaleMarginalLikl_mat_ls, rowSums))


  if(V == 1){
    LL_GVmat2 = LL_GVmat #adjust for prior
    postGV = matrix(1,nrow = NG, ncol = V)
    row.names(postGV) = row.names(std.data)
  }else{
    LL_GVmat2 = t(apply(LL_GVmat,1,function(x) log(pV)+x)) #adjust for prior
    baseLL = sapply(1:V, function(v) {
      K_v = K[v]
      pvK = pVK[[v]]
      mu_GK_v = mu_GK_V[[v]]
      sigma_Gv = sigma_GV[,v]

      if(K_v>1){
        idx_pair = combn(K_v,2)
        diffmu2G = apply(idx_pair,2,function(xx) {
          a = xx[1]
          b = xx[2]
          (pvK[a] * pvK[b])*(mu_GK_v[,a]-mu_GK_v[,b])^2
        })
        scale_p = sum(apply(idx_pair,2,function(xx) {
          a = xx[1]
          b = xx[2]
          (pvK[a] * pvK[b])
        }))
        rowSums(diffmu2G)/(scale_p*2*sigma_Gv^2) + log(sigma_Gv*sqrt(2*pi))
      }else{
        rep(0,length(sigma_Gv))
      }
    })
    LL_GVmat3 = LL_GVmat2 + N*baseLL #adjust for relative likelihood
    postGV = t(apply(LL_GVmat3,1, function(x) {
      exp(x-max(x))/sum(exp(x-max(x)))
    }))
  }
  logScaleConditionalLikl_ls = ConditionalLL_ls_func(N, NG, V, K, std.data,
                                                     pVK, mu_GK_V, sigma_GV)
  w_NK_V = lapply(1:V, function(v){
    w_Gv = postGV[,v]
    logScaleConditionalLikl_v = logScaleConditionalLikl_ls[[v]]
    wlogScaleConditionalLikl_v = sapply(logScaleConditionalLikl_v,function(xx){
      apply(xx,2,function(yy) sum(yy*w_Gv))
    })
    w_NK_v = t(sapply(1:nrow(wlogScaleConditionalLikl_v), function(i) {
      logScale = wlogScaleConditionalLikl_v[i,]
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
    postGV[nonclust_idx,] = 0
    clust_idx = which(maxR2 >= R2_cutoff)
    predGV_clust = apply(postGV[clust_idx,], 1, which.max)

  }else{
    R2_V = as.matrix(1-SSE_GV/SSE_G0)
    maxR2 = R2_V[,1]
    nonclust_idx = which(maxR2 < R2_cutoff)
    postGV[nonclust_idx,] = 0
    clust_idx = which(maxR2 >= R2_cutoff)
    predGV_clust = rep(1,length(clust_idx))

  }

  tb = table(predGV_clust)
  null.v = union(as.numeric(names(tb)[tb<4]),which(apply(postGV,2,function(x) sum(x > 1e-10) == 0)))
  if(length(null.v) != 0){
    V = V-length(null.v)
    if(V == 0){
      print("Clusterable views are reduced to 0, consider a smaller R2 cutoff")
      return(list(V=V, K=K, pV=pV, pVK=pVK, mu_GK_V=mu_GK_V,  sigma_GV =sigma_GV,
                  Kpath=Kpath,postGV = w_GV_star, postNK_V = w_NK_V,
                  R2_cutoff = R2_cutoff,
                  avgLLik_selected_soft_sepV=NA,avgR2_selected_soft_sepV = NA,
                  criterion = criterion))
    }
    K = K[-null.v]
    postGV = postGV[,-(null.v)]
    idx = apply(postGV,1,function(x){!all(x == 0)})
    if(ncol(postGV) == 1){
      postGV[idx,] = apply(postGV[idx,,drop = F], 1, function(xx) xx/sum(xx))
    }else{
      postGV[idx,] = t(apply(postGV[idx,], 1, function(xx) xx/sum(xx)))
    }
    w_NK_V = w_NK_V[-null.v]
    mu_GK_V = mu_GK_V[-null.v]
    sigma_GV = sigma_GV[,-null.v,drop = F]
    pV = pV[-(null.v)]
    pV = pV/(sum(pV))
    pVK = pVK[-null.v]
    LL_GVmat = LL_GVmat[,-null.v,drop = F]
    R2_V = R2_V[,-null.v,drop = F]

    V = length(mu_GK_V)
    print(paste0("Adjust V to ", V,"!"))
  }
  if(V != 1){
    postGV_hard = t(apply(postGV, 1, function(xx) {
      out = rep(0,length(xx))
      out[which.max(xx)] = 1
      return(out)
    }))
  }else{
    postGV_hard = postGV
  }
  Lik_Gmat = rowSums(postGV*LL_GVmat)
  avgLLik_selected_soft_sepV = mean(sapply(1:V, function(v){
    xx = which(apply(postGV, 1, which.max) == v)
    mean(Lik_Gmat[intersect(xx,clust_idx)])
  }))

  avgR2_selected_soft_sepV = mean(sapply(1:V, function(v){
    xx = which(apply(postGV, 1, which.max) == v)
    mean(rowSums(postGV * R2_V)[intersect(xx,clust_idx)])
  }))
  return(list(V=V, K=K, pV=pV, pVK=pVK, mu_GK_V=mu_GK_V,sigma_GV =sigma_GV,
              Kpath=Kpath,postGV = postGV, postNK_V = w_NK_V,
              R2_cutoff = R2_cutoff,
              avgLLik_selected_soft_sepV=avgLLik_selected_soft_sepV,avgR2_selected_soft_sepV = avgR2_selected_soft_sepV,
              criterion = criterion))
}

