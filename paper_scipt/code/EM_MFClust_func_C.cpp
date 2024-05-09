#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}


// [[Rcpp::export]]
List mat_to_list(NumericMatrix m, int V, NumericVector K) {
  List l(V);
  for (int v = 0; v < V; ++v) {
    int start_idx;
    int end_idx;
    if(v == 0){
      start_idx = 0;
      end_idx = (K[v]-1);
    }else{
      IntegerVector r1 = seq(0,v-1);
      IntegerVector r2 = seq(0,v);

      start_idx = sum(as<NumericVector>(K[r1]));
      end_idx = sum(as<NumericVector>(K[r2]))-1;
    }
    l[v] = m(_,Range(start_idx,end_idx));
  }
  return l;
}

// [[Rcpp::export]]
NumericMatrix sigma_GKV_func(List mu_GK_V_new, NumericMatrix std_data, List w_NK_V,
                             int N, int NG, int V){
  NumericMatrix sigma_GKV_new(NG,V);
  for(int v = 0; v < V; ++v){
    NumericMatrix mu_GK_v = mu_GK_V_new[v];
    NumericMatrix w_NK_v = w_NK_V[v];

    for(int g = 0; g < NG; ++g){
      double asum_i = 0;
      for (int k = 0; k < mu_GK_v.ncol(); ++k) {
        NumericVector centerY = std_data(g,_) - mu_GK_v(g,k);
        NumericVector w_k = w_NK_v(_,k);
        asum_i += sum(w_k*centerY*centerY);
      }
      sigma_GKV_new(g,v) = pow(asum_i/N,0.5);
    }
  }
  return sigma_GKV_new;
}
