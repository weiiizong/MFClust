# MFClust
This packge implemetn the model-based Multi-Facet Clustering (MFClust) algorithm to simultaneously identify multiple facets of potential clustering solutions/facets in high dimensional omics data. Two human postmortem brain tissue datasets after preprocessing that have been used in the manuscript are included: 1) miacroarray data from two brain regions, Brodmann’s area 11 (BA11) and BA47, 2) RNA-Seq data from the nucleus accumbens. Two main functions are included: *MFClust_init* applies the tight clustering initialization algorithm to generate reasonable initials and *MFClust_init* implements the model based multi-F=facet Clustering (MFClust) algorithm given initials.

## Install
MFClust package files are in MFClust/ folder, You can install by copying and paste the following code into R
```
devtools::install_github("weiiizong/MFClust/MFClust")
```
Alternatively, download the tar.gz zipped file, unzip it and install using devtools::install(). Make sure R packages Rcpp, mclust, tightClust, cluster, MASS have all been properly imported.
