# MFClust
This packge implemetn the model-based Multi-Facet Clustering (MFClust) algorithm to simultaneously identify multiple facets of potential clustering solutions/facets in high dimensional omics data. Two human postmortem brain tissue datasets after preprocessing that have been used in the manuscript are included: 1) miacroarray data from two brain regions, Brodmann’s area 11 (BA11) and BA47, 2) RNA-Seq data from the nucleus accumbens. Two main functions are included: *MFClust_init* applies the tight clustering initialization algorithm to generate reasonable initials and *MFClust_init* implements the model based multi-F=facet Clustering (MFClust) algorithm given initials.

## Installation
MFClust package files are in MFClust/ folder, You can install by copying and paste the following code into R
```
devtools::install_github("weiiizong/MFClust/MFClust")
```
Alternatively, download the tar.gz zipped file, unzip it and install using devtools::install(). Make sure R packages Rcpp, mclust, tightClust, cluster, MASS have all been properly imported.

## Initialization
```
data(BA11_BA47_NG2000_dat)
init.ls = MFClust_init(std.data=BA11_BA47_NG2000_dat, V=5, G.K=20, initial.kmin=100,
R2_cutoff = 0.26, seed = 9)
```

## Application of algorithm
```
#extract parameters from init.ls
V = length(init.ls)
pV=rep(1/V,V)
K=sapply(init.ls[1:V], function(xx) ncol(xx$mean.mat))
pVK=lapply(1:V, function(v) rep(1/K[v],K[v]))
mu_GK_V=lapply(init.ls[1:V], function(xx) {
mat = xx$mean.mat
row.names(mat) = row.names(BA11_BA47_NG2000_dat)
return(mat)
})
sigma_GV=sapply(init.ls[1:V], function(xx) {
  mat = apply(xx$sigma.mat,1,function(x){sqrt(mean(x^2))})
  return(mat)
})
row.names(sigma_GV) = row.names(BA11_BA47_NG2000_dat)
res = MFClust(V, K, pV, pVK, mu_GK_V, sigma_GV, BA11_BA47_NG2000_dat, R2_cutoff = 0.26)
```
