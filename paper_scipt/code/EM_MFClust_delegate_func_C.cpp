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
IntegerVector which_maxCpp(NumericVector v) {
  
  double current_max = v[0];
  int n = v.size();
  std::vector< int > res;
  res.push_back( 0 );
  int i;
  
  for( i = 1; i < n; ++i) {
    double x = v[i];
    if( x > current_max ) {
      res.clear();
      current_max = x;
      res.push_back( i );
    } else if ( x == current_max ) {
      res.push_back( i );
    }
  }
  Rcpp::IntegerVector iv( res.begin(), res.end() );
  return iv;
}

// [[Rcpp::export]]
NumericVector logLik_vG(int N, int NG, int K_v, double pv, NumericMatrix std_data, //f_gv
                                NumericMatrix mu_GK, NumericVector sigma_Gv, NumericVector pK) {
  NumericVector logLik(NG);
  
  for(int g = 0; g < NG; ++g){
    NumericVector logLik_vgN(N);
    for(int n=0; n < N; ++n){
      double y = std_data(g,n);
      NumericVector xx(K_v);
      for(int k = 0; k < K_v; ++k){
        xx[k] = R::dnorm(y, mu_GK(g,k), sigma_Gv[g], false);
      }
      NumericVector xx2 = pK*xx;
      NumericVector t_gn = xx2/sum(xx2);
      logLik_vgN[n] = sum(t_gn*(log(xx2)));
    }
    logLik[g] = log(pv) + sum(logLik_vgN);
  }
  return logLik;
}

// [[Rcpp::export]]
NumericVector logLik_vG_hard(int N, int NG, int K_v, double pv, NumericMatrix std_data, //f_gv
                        NumericMatrix mu_GK, NumericVector sigma_Gv, NumericVector pK) {
  NumericVector logLik(NG);
  
  for(int g = 0; g < NG; ++g){
    NumericVector logLik_vgN(N);
    for(int n=0; n < N; ++n){
      double y = std_data(g,n);
      NumericVector xx(K_v);
      for(int k = 0; k < K_v; ++k){
        xx[k] = R::dnorm(y, mu_GK(g,k), sigma_Gv[g], false);
      }
      NumericVector xx2 = pK*xx;
      int idx = which_maxCpp(xx2)[0];
      logLik_vgN[n] = log(xx2)[idx];
    }
    logLik[g] = log(pv) + sum(logLik_vgN);
  }
  return logLik;
}

// [[Rcpp::export]]
NumericVector logLik_vG_hard2(int N, int NG, int K_v, double pv, NumericMatrix std_data, //f_gv
                             NumericMatrix mu_GK, NumericVector sigma_Gv) {
  NumericVector logLik(NG);
  
  for(int g = 0; g < NG; ++g){
    NumericVector logLik_vgN(N);
    for(int n=0; n < N; ++n){
      double y = std_data(g,n);
      NumericVector xx(K_v);
      for(int k = 0; k < K_v; ++k){
        xx[k] = R::dnorm(y, mu_GK(g,k), sigma_Gv[g], false);
      }
      NumericVector xx2 = xx;
      int idx = which_maxCpp(xx2)[0];
      logLik_vgN[n] = log(xx2)[idx];
    }
    logLik[g] = log(pv) + sum(logLik_vgN);
  }
  return logLik;
}

// [[Rcpp::export]]
List ConditionalLL_ls_func(int N, int NG, int V, NumericVector K, NumericMatrix std_data,
                           List mu_GK_V, NumericMatrix sigma_GV) {
  List ConditionalLL(V);
  for(int v = 0; v < V; ++v){
    int K_v = K[v];
    NumericMatrix mu_GK =  as<NumericMatrix>(mu_GK_V[v]);

    List ConditionalLikl_v(K_v);
    for(int k = 0; k < K_v; ++k){
      NumericMatrix Likl_k(NG,N);
      for(int g = 0; g < NG; ++g){
        NumericVector y = std_data(g,_);
        // Rcout << "y = "<<y<<"\n";
        
        Likl_k(g,_) = Rcpp::dnorm(y, mu_GK(g,k), sigma_GV(g,v), true);
        // Rcout << "mu_GK(g,k) = "<<mu_GK(g,k)<<"\n";
        // Rcout << "sigma_GK(g,k) = "<<sigma_GK(g,k)<<"\n";
        // 
        // Rcout << "Likl_k(g,_) = "<<NumericVector{Likl_k(g,_)}<<"\n";
        
      }
      ConditionalLikl_v[k] = Likl_k;
    }
    ConditionalLL[v] = ConditionalLikl_v;
  }
  return ConditionalLL;
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
