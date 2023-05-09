run.func_V = function(l){
  initial.idx = scenario.idx[l,1]
  R2.idx = scenario.idx[l,2]
  data.idx = scenario.idx[l,2]
  
  initial.kmin = initial.list[initial.idx,"kmin"]
  initial.seednum = initial.list[initial.idx,"seednum"]
  
  R2_cutoff = R2.list[R2.idx]
  std.data = data.list[[data.idx]]
  mu0 = mu0.list[[data.idx]]
  V = V.list[[scenario.idx[l,3]]]
  
  #initialization =======================
  data = std.data
  set.seed(initial.seednum)
  
  out = tryCatch({
    G.K = 20
    tClust = tight.clust(data, target = G.K, k.min = initial.kmin)
    G.K0 = as.numeric(setdiff(names(table(tClust$cluster))[table(tClust$cluster) > 3],"-1"))
    if(length(G.K0) == 0){
      init.ls = list()
      aout = list(res = NA,res.tb = NA,metric.tb = NA)
    }else if(length(G.K0) < 4){ #not enough for module clustering
      init.ls = lapply(G.K0, function(v){
        sub.data = data[which(tClust$cluster == v),]
        
        #refit kmeans using combined genes
        res.S4 = K.Clust(
          t(sub.data),
          method = "S4",
          Kmin = 2,
          Kmax = 6,
          trim.S4 = 0.05,
          cutoff = 0,
          n.resample = 50
        )
        res.Kmeans = kmeans(t(sub.data),centers=res.S4[[1]],nstart = 100)
        clust.lb_m = res.Kmeans$cluster
        
        
        clust.fit.ls = lapply(1:res.S4[[1]], function(k){
          clust.expr = std.data[,clust.lb_m == k,drop = F]
          clust.mean = apply(clust.expr,1, function(m) {
            fitdistr(m, "normal")$estimate["mean"]
          })
          clust.sigma = apply(clust.expr,1, function(m) {
            fitdistr(m, "normal")$estimate["sd"]
          })
          clust.lb = as.numeric(clust.lb_m == k)
          return(list(clust.mean = clust.mean, clust.sigma = clust.sigma, clust.lb = clust.lb))
        })
        mean.mat = sapply(1:length(clust.fit.ls), function(x) clust.fit.ls[[x]]$clust.mean)
        sigma.mat = sapply(1:length(clust.fit.ls), function(x) clust.fit.ls[[x]]$clust.sigma)
        lb.mat = sapply(1:length(clust.fit.ls), function(x) clust.fit.ls[[x]]$clust.lb)
        pVk = apply(lb.mat,2,sum)/sum(apply(lb.mat,2,sum))
        
        list(sub.data = sub.data, mean.mat = mean.mat, sigma.mat = sigma.mat,cl.lb = clust.lb_m,
             lb.mat = lb.mat, gene.lb = which(tClust$cluster == v), pVk = pVk)
      })
      
    }else{
      gene.module.res = lapply(G.K0,function(tc){
        print(tc)
        sub.data = data[which(tClust$cluster == tc),]
        res.S4 = K.Clust(
          t(sub.data),
          method = "S4",
          Kmin = 2,
          Kmax = 6,
          trim.S4 = 0.05,
          cutoff = 0,
          n.resample = 50
        )
        res.Kmeans = kmeans(t(sub.data),centers=res.S4[[1]],nstart = 100)
        clust.lb_m = res.Kmeans$cluster
        
        #Evaluate R2 for genes in each module
        aclust.fit.ls = lapply(1:res.S4[[1]], function(k){
          clust.expr = sub.data[,clust.lb_m == k,drop = F]
          clust.mean = apply(clust.expr,1, function(m) {
            fitdistr(m, "normal")$estimate["mean"]
          })
          clust.lb = as.numeric(clust.lb_m == k)
          return(list(clust.mean = clust.mean, clust.lb = clust.lb))
        })
        amean.mat = sapply(1:length(aclust.fit.ls), function(x) aclust.fit.ls[[x]]$clust.mean)
        alb.mat = sapply(1:length(aclust.fit.ls), function(x) aclust.fit.ls[[x]]$clust.lb)
        
        sub.data.pred = amean.mat %*% t(alb.mat)
        SSE_G = apply((sub.data.pred-sub.data)^2, 1, sum)
        
        pred_GN0 = matrix(rep(apply(sub.data, 1, mean),ncol(sub.data)),nrow = nrow(sub.data),ncol = ncol(sub.data))
        SSE_G0 = apply((pred_GN0-sub.data)^2, 1, sum)
        R2 = 1-SSE_G/SSE_G0
        
        if(sum(R2< R2_cutoff) != 0){
          idx = which(tClust$cluster == tc)[-which(R2< R2_cutoff)]
          
          if(length(idx) <= 1){
            return(NULL)
          }
          #refit clusters
          sub.data = data[idx,]
          res.S4 = K.Clust(
            t(sub.data),
            method = "S4",
            Kmin = 2,
            Kmax = 6,
            trim.S4 = 0.05,
            cutoff = 0,
            n.resample = 50
          )
          res.Kmeans = kmeans(t(sub.data),centers=res.S4[[1]],nstart = 100)
          clust.lb_m = res.Kmeans$cluster
          
        }else{
          idx = which(tClust$cluster == tc)
        }
        
        
        return(list(clust.lb=clust.lb_m, idx=idx))
      })
      gene.module.res = gene.module.res[!sapply(gene.module.res,is.null)]
      clust.lb.ls = lapply(gene.module.res, "[[","clust.lb")
      gene.lb.ls = lapply(gene.module.res, "[[","idx")
      
      pairwiseARImat = sapply(1:length(clust.lb.ls), function(x1){
        rm1 = setdiff(1:length(clust.lb.ls),x1)
        sapply(1:length(clust.lb.ls), function(x2){
          adjustedRandIndex(clust.lb.ls[[x1]],clust.lb.ls[[x2]])
        })
      })
      distARImat = 1-pairwiseARImat
      
      clust.genes = pam(distARImat, V, diss = TRUE)
      
      clust.genes.lb = clust.genes$clustering
      
      init.ls = lapply(1:V, function(v){
        sub.data = data[unlist(gene.lb.ls[clust.genes.lb == v]),]
        
        #refit kmeans using combined genes
        res.S4 = K.Clust(
          t(sub.data),
          method = "S4",
          Kmin = 2,
          Kmax = 6,
          trim.S4 = 0.05,
          cutoff = 0,
          n.resample = 50
        )
        res.Kmeans = kmeans(t(sub.data),centers=res.S4[[1]],nstart = 100)
        clust.lb_m = res.Kmeans$cluster
        
        
        clust.fit.ls = lapply(1:res.S4[[1]], function(k){
          clust.expr = std.data[,clust.lb_m == k,drop = F]
          clust.mean = apply(clust.expr,1, function(m) {
            fitdistr(m, "normal")$estimate["mean"]
          })
          clust.sigma = apply(clust.expr,1, function(m) {
            fitdistr(m, "normal")$estimate["sd"]
          })
          clust.lb = as.numeric(clust.lb_m == k)
          return(list(clust.mean = clust.mean, clust.sigma = clust.sigma, clust.lb = clust.lb))
        })
        mean.mat = sapply(1:length(clust.fit.ls), function(x) clust.fit.ls[[x]]$clust.mean)
        sigma.mat = sapply(1:length(clust.fit.ls), function(x) clust.fit.ls[[x]]$clust.sigma)
        lb.mat = sapply(1:length(clust.fit.ls), function(x) clust.fit.ls[[x]]$clust.lb)
        pVk = apply(lb.mat,2,sum)/sum(apply(lb.mat,2,sum))
        
        list(sub.data = sub.data, mean.mat = mean.mat, sigma.mat = sigma.mat,cl.lb = clust.lb_m,
             lb.mat = lb.mat, gene.lb = unlist(gene.lb.ls[clust.genes.lb == v]), pVk = pVk)
      })
    }
    
    #run algorithm =======================
    if(length(init.ls) != 0){
      V = length(init.ls)
      pV=rep(1/V,V)
      K=sapply(init.ls[1:V], function(xx) ncol(xx$mean.mat))
      pVK=lapply(1:V, function(v) rep(1/K[v],K[v]))
      mu_GK_V=lapply(init.ls[1:V], function(xx) {
        mat = xx$mean.mat
        row.names(mat) = row.names(std.data)
        return(mat)
      })
      sigma_GV=sapply(init.ls[1:V], function(xx) {
        mat = apply(xx$sigma.mat,1,function(x){sqrt(mean(x^2))})
        return(mat)
      })
      row.names(sigma_GV) = row.names(std.data)
      
      res = EM_multiView(V, K, pV, pVK, mu_GK_V, sigma_GV, std.data,R2_cutoff = R2_cutoff,
                         quite=T, updateK = T, updateK_thin = 1, maxr = 200)
      if(res$V < 2){
        res.tb = c(mu0 = mu0, R2_cutoff=R2_cutoff, orig_V = V, initial.kmin=initial.kmin, initial.seednum=initial.seednum,
                   max_pairARI=NA,avg_pairARI=NA,V=res$V,NGclust=0,
                   avgR2_selected_soft = res$avgR2_selected_soft, avgR2_selected_hardView = res$avgR2_selected_hardView,avgR2_selected_hardViewClust=res$avgR2_selected_hardViewClust,
                   minR2_selected_soft = res$minR2_selected_soft, minR2_selected_hardView = res$minR2_selected_hardView, minR2_selected_hardViewClust = res$minR2_selected_hardViewClust,
                   avgR2_selected_soft_sepV = res$avgR2_selected_soft_sepV, avgR2_selected_hardView_sepV = res$avgR2_selected_hardView_sepV,avgR2_selected_hardViewClust_sepV=res$avgR2_selected_hardViewClust_sepV,
                   minR2_selected_soft_sepV = res$minR2_selected_soft_sepV, minR2_selected_hardView_sepV = res$minR2_selected_hardView_sepV, minR2_selected_hardViewClust_sepV = res$minR2_selected_hardViewClust_sepV,
                   initial.idx=initial.idx,R2.idx=R2.idx,data.idx=data.idx,l=l)
        metric.tb = NA
      }else{
        cl.ls = lapply(1:res$V, function(v){
          apply(res$w_NK_V[[v]],1,which.max)
        })
        
        Vpairs = combn(res$V,2)
        sample_Clust_ARI = sapply(1:ncol(Vpairs), function(z) {
          mclust::adjustedRandIndex(factor(cl.ls[[Vpairs[1,z]]]),factor(cl.ls[[Vpairs[2,z]]]))
        })
        max_pairARI = max(sample_Clust_ARI)
        avg_pairARI = mean(sample_Clust_ARI)
        
        pred_Glb = apply(res$postGV,1,function(x) {
          if(all(x == 0)){
            0
          }else{
            which.max(x)
          }
        })
        res.tb = c(mu0 = mu0, R2_cutoff=R2_cutoff, orig_V = V, initial.kmin=initial.kmin, initial.seednum=initial.seednum,
                   max_pairARI=max_pairARI,avg_pairARI=avg_pairARI,V=res$V,NGclust=sum(pred_Glb != 0),
                   avgR2_selected_soft = res$avgR2_selected_soft, avgR2_selected_hardView = res$avgR2_selected_hardView,avgR2_selected_hardViewClust=res$avgR2_selected_hardViewClust,
                   minR2_selected_soft = res$minR2_selected_soft, minR2_selected_hardView = res$minR2_selected_hardView, minR2_selected_hardViewClust = res$minR2_selected_hardViewClust,
                   avgR2_selected_soft_sepV = res$avgR2_selected_soft_sepV, avgR2_selected_hardView_sepV = res$avgR2_selected_hardView_sepV,avgR2_selected_hardViewClust_sepV=res$avgR2_selected_hardViewClust_sepV,
                   minR2_selected_soft_sepV = res$minR2_selected_soft_sepV, minR2_selected_hardView_sepV = res$minR2_selected_hardView_sepV, minR2_selected_hardViewClust_sepV = res$minR2_selected_hardViewClust_sepV,
                   initial.idx=initial.idx,R2.idx=R2.idx,data.idx=data.idx,l=l)
        metric.tb = t(sapply(1:res$V, function(v){
          cl = apply(res$w_NK_V[[v]],1,which.max)
          ari_res = sapply(1:length(clust.lb0), function(xx){
            adjustedRandIndex(clust.lb0[[xx]],cl)
          })
          Cl_res = c(ARI = max(ari_res), which = which.max(ari_res))
          orig.idx = which.max(ari_res)
          
          binary_lb = factor(as.numeric(pred_Glb == v))
          true.gene.lb = rep(0, nrow(std.data))
          true.gene.lb[gene.lb0[[orig.idx]]] = 1
          true.gene.lb = factor(true.gene.lb)
          mat = confusionMatrix(binary_lb,true.gene.lb,positive = "1")
          G_Sensitivity = mat$byClass["Sensitivity"]
          G_Specificity = mat$byClass["Specificity"]
          
          return(c(v = v, orig.idx = orig.idx, Cl_ARI = max(ari_res), NG = sum(pred_Glb == v),
                   G_Sensitivity, G_Specificity))
        }))
      }
      aout = list(res = res,res.tb = res.tb,metric.tb = metric.tb)
    }
    aout
  },error = function(e) {
    return(list(e = e, l=l))
  })
  return(out)
}
