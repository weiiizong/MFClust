S4_Kmedoids<-function(x, K.try = 2:10,n.resample,trim=0.05,cut.off=0.8,n.prop = 0.7){
  n = nrow(x)
  if(n >= 9){
    sub.n = round(n*n.prop)
  }else{
    sub.n = n-1
  }
  result_k<-function(x,K.try){
    result = sapply(K.try, function(K){
      print(K)
      #o.cluster <- kmeans(x, centers = K, nstart = 100)$cluster
      result.subsampling = lapply(1:n.resample, function(b){
        index.subsample = sample(1:nrow(x), sub.n, replace = F)
        #index.subsample = unlist(sapply(1:K,function(i) sample(which(o.cluster == i),
        #                         round(sum(o.cluster == i)*0.7), replace = F)))
        xb = x[index.subsample,index.subsample]
        km.out <- pam(xb, k = K, diss  = TRUE)
        # run sparse k-means
        Cb = km.out$cluster
        
        group.init = rep(NA, n)
        group.init[index.subsample] = Cb
        consensus.matrix = sapply(1:n, function(i){
          if(i %in% index.subsample){
            as.integer(group.init[i] == group.init)##1==NA is NA
          } else rep(NA, n)
        })
        #consensus.matrix.upper = consensus.matrix[upper.tri(consensus.matrix)]
        return(list(clustering = consensus.matrix))
      })
      
      cluster = sapply(result.subsampling, function(x) x$clustering, simplify = "array")
      
      r.mtx  =  apply(cluster, c(1,2), function(x) mean(x != 0, na.rm = T))###the average concensus matrix of subsample
      
      o.cluster <- pam(x, k = K, diss = TRUE)$cluster####kmeans to the whole sample
      
      o.mtx = sapply(1:n, function(i){
        as.integer(o.cluster[i] == o.cluster)
      })###the concensus matrix of whole sample
      
      r.mtx[which(o.mtx == 0)] <- (1 - r.mtx[which(o.mtx == 0)])
      nr<-nrow(x)
      rm1.calc = function(index){
        rm1 = vector("numeric")
        r.mtx1 = r.mtx[index,index]
        o.mtx1 = o.mtx[index,index]
        for(i in 1:length(index)){
          if(all(c(1,0)%in%o.mtx1[i,])){
            rm1[i] = (mean(r.mtx1[i, which(o.mtx1[i,] == 1)],na.rm = T) +
                        mean(r.mtx1[i, which(o.mtx1[i,] == 0)],na.rm = T))-1
          } else if(1%in%o.mtx1[i,]){
            rm1[i] = mean(r.mtx1[i, which(o.mtx1[i,] == 1)],na.rm = T)
          } else {
            rm1[i] = mean(r.mtx1[i, which(o.mtx1[i,] == 0)],na.rm = T)
          }
        }
        return(rm1)
      }
      # #res_score<-rev(sort(rm1.calc(1:nr)))
      for(i in 1:(nr-1)){
        index.order<-order(rm1.calc(i:nr))
        r.mtx[i:nr,i:nr] = r.mtx[i:nr,i:nr][index.order,
                                            index.order]
        #index[i:nr]<-index[i:nr][index.order]
        o.mtx[i:nr,i:nr] = o.mtx[i:nr,i:nr][index.order,
                                            index.order]
      }
      Num<-round(trim*nr)
      r.mtx<-r.mtx[(Num+1):nr,(Num+1):nr]
      o.mtx<-o.mtx[(Num+1):nr,(Num+1):nr]
      stat.calc1 = function(i){
        nr1<-dim(o.mtx)[1]
        if(all(c(1,0)%in%o.mtx[i,1:nr1])){
          (mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 0)]) +
             mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 1)])) - 1
        }
        else if(1%in%o.mtx[i,1:nr1]){
          mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 1)])
        }
        else {
          mean(r.mtx[i,1:nr1][which(o.mtx[i,1:nr1] == 0)])
        }
      }
      
      #per<-1-trim
      stat = rev(sapply(1:(dim(o.mtx)[1]), stat.calc1))
      stat<-sort(stat,decreasing = T)
      # plot(stat,ylim = c(0,1))
      # #title(paste("plot of S4.Iter","_K=",K,sep=""))
      # #barplot(stat,width = 0.832)
      # title(paste("k=",K,sep=""))
      # abline(v=n*0.95,col="red",lwd=3.5)
      U.min = mean(stat)
      return(U.min)
    })
  }
  result_data<-result_k(x,K.try)
  
  
  #if(specificity==TRUE){
  #  result<-apply(result_data,2,mean)
  #}else{
  #  result<-result_data[1,]
  #}
  #if((res.null==1)&(result_data==1)){
  #  maxk<-K.try[max(which(result_data == max(result_data)))]
  #}else{
  if(length(which(is.na(result_data)))!=0){
    result_data[which(is.na(result_data))]<-0
  }
  if(max(result_data)<cut.off){
    maxk<-1
  }else{
    maxk<-K.try[max(which(result_data == max(result_data)))]
  }
  
  #}
  res<-list(maxk,score=result_data)
  return(res)
}
