#' Initialization of the model-based multi-facet clustering algorithm
#'
#' @param std.data a standardized data matrix. Rows are the genes/features and columns are
#' samples. Each row is standardized to mean zero and standard deviation one.
#' @param V the number of facets/clustering solutions.
#' @param G.K the total number of clusters that the user aims to find in tight clustering algorithm.
#' @param initial.kmin the starting point of k0 in tight clustering algorithm.
#' @param R2_cutoff genes/features whose fitted R2 is smaller than this value will be labeled
#' as non-clusterable genes and do not belong to any facet.
#' @param seed random seed number. Default is NA.
#'
#' @return A list of initial parameters in each facet
#' @export MFClust_init
#'
#' @examples
#' \dontrun{
#' data(BA11_BA47_NG2000)
#' init.ls = MFClust_init(std.data=BA11_BA47_NG2000, V=5, G.K=20, initial.kmin=100,
#' R2_cutoff = 0.26, seed = 9)
#' }
MFClust_init = function(std.data, V, G.K=20, initial.kmin=100, R2_cutoff, seed = NA){
  if(!is.na(seed)){
    set.seed(seed)
  }
  tClust = tight.clust(std.data, target = G.K, k.min = initial.kmin)
  G.K0 = as.numeric(setdiff(names(table(tClust$cluster))[table(tClust$cluster) > 3],"-1"))
  if(length(G.K0) == 0){
    init.ls = list()
    aout = list(res = NA,res.tb = NA,metric.tb = NA)
  }else if(length(G.K0) < V){ #not enough for module clustering
    init.ls = lapply(G.K0, function(v){
      sub.data = std.data[which(tClust$cluster == v),]

      #refit kmeans using combined genes
      res.S4 = K.Clust(
        t(sub.data),
        method = "S4",
        Kmin = 2,
        Kmax = 5,
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
      sub.data = std.data[which(tClust$cluster == tc),]
      res.S4 = K.Clust(
        t(sub.data),
        method = "S4",
        Kmin = 2,
        Kmax = 5,
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
        sub.data = std.data[idx,]
        res.S4 = K.Clust(
          t(sub.data),
          method = "S4",
          Kmin = 2,
          Kmax = 5,
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
      sub.data = std.data[unlist(gene.lb.ls[clust.genes.lb == v]),]

      #refit kmeans using combined genes
      res.S4 = K.Clust(
        t(sub.data),
        method = "S4",
        Kmin = 2,
        Kmax = 5,
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
}
