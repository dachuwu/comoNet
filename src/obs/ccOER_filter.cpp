#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

dcube proc_OERmat(sp_mat A){
  
  sp_mat OER(A.n_rows, A.n_cols);
  sp_mat sdlOER(A.n_rows, A.n_cols);
  int nn = A.n_rows;
  double S = accu(A);
  sp_mat s_i = sum(A, 1);
  sp_mat s_j = sum(A, 0);
  
  for(int i=0; i < nn; ++i){
    for(int j=0; j < nn; ++j){
      double a = A(i,j)/S;
      if(a==0) continue;
      double b = s_i(i,0)/S - a;
      double c = s_j(0,j)/S - a;
      OER(i,j) = a/((a+b)*(a+c));
      sdlOER(i,j) = (b*c*(1-a) + (1-b-c)*pow(a,2) - pow(a,3))/(a*(a+b)*(a+c)*S);
    }
  }
  sdlOER = sqrt(sdlOER);
  
  dcube out(A.n_rows, A.n_cols, 2);
  out.slice(0) = OER;
  out.slice(1) = sdlOER;
  return out;
}

List proc_OERvec(sp_mat A, IntegerVector eval_i, IntegerVector eval_j){
  
  int pp = eval_i.length();
  NumericVector OER ;
  NumericVector sdlOER ;
  double S = accu(A);
  sp_mat s_i = sum(A, 1);
  sp_mat s_j = sum(A, 0);
  for(int p=0; p < pp; ++p){
    int i = eval_i[p]; 
    int j = eval_j[p];
    
    if(A(i,j)==0){
      OER.push_back(0);
      sdlOER.push_back(0);
    } else {
      double a = A(i,j)/S;
      double b = s_i(i,0)/S - a;
      double c = s_j(0,j)/S - a;
      OER.push_back( a/((a+b)*(a+c)));
      sdlOER.push_back((b*c*(1-a) + (1-b-c)*pow(a,2) - pow(a,3))/(a*(a+b)*(a+c)*S));
    }
    
  }
  sdlOER = sqrt(sdlOER);
  return List::create(_["OER"]= OER,  _["sdlOER"]= sdlOER );
  
}

void uvec_push(uvec & v, unsigned int value) {
  uvec av(1);
  av(0) = value;
  v.insert_rows(v.n_rows, av);
}


List BDC_cs(sp_mat Asp){
  
  dmat A(Asp);
  int nn = A.n_rows;
  List out;
  for(int i=0; i < nn; ++i){
    std::vector<int> nei;
    for(int j=0; j < nn; ++j){
      if(A(i,j) >0 ) nei.push_back(j);
    }
    
    if(nei.size()>1){
      uvec nei2 = conv_to<uvec>::from(nei);
      uvec rk;
      uvec ri;
      uvec rj;
      for(std::vector<int>::iterator it = nei.begin() ; it != nei.end(); ++it){
        uvec sr(1); sr[0]= *it;
        uvec lg_nei = find( A.submat(sr, nei2) > 0);
        if(lg_nei.n_rows > 0){
          uvec k_sub(lg_nei.n_rows); k_sub.fill(i);
          uvec i_sub(lg_nei.n_rows); i_sub.fill(*it);
          rk = join_cols(rk, k_sub );
          ri = join_cols(ri, i_sub );
          rj = join_cols(rj, nei2(lg_nei));
        }
      }
      
      if(rk.n_elem > 0){
        uvec_push(rk, 0);
        uvec_push(ri, 0);
        uvec_push(rj, 0);
        out.push_back(DataFrame::create(_["k"]=rk, _["i"]=ri, _["j"]=rj));}
    }
  }
  return out;
}


// [[Rcpp::export]]
List filter_ccOER(arma::sp_mat W, List U, double sig_lv){
  
  double zs =  R::qnorm(sig_lv/2, 0.0, 1.0, 1, 0);   /////// -1.96
  
  dcube est0 = proc_OERmat(W);
  sp_mat lgW_0( est0.slice(0) % exp(est0.slice(1)*zs) );
  
  lgW_0.transform([](double val){ if(val>1.0){ return 1.0; }else{ return 0.0;} } );
  lgW_0.clean(0.1);
  sp_mat lgW_1(lgW_0);
  List kij = BDC_cs(W % lgW_0);
  
  IntegerVector rec_k;
  IntegerVector rec_i;
  IntegerVector rec_j;
  
  for(int it = 0 ; it < kij.size() ; ++it){
    
    DataFrame dd = kij[it];
    
    if(dd.nrows() >= 2 ){
      IntegerVector dd_k = dd["k"];
      IntegerVector dd_i = dd["i"];
      IntegerVector dd_j = dd["j"];
      
      sp_mat W_k = U[dd_k[1]];
      
      List res0 = proc_OERvec(W - W_k, dd_i, dd_j);
      NumericVector res0_oer = res0["OER"];
      NumericVector res0_sdloer = res0["sdlOER"];
      NumericVector cc0_lb = res0_oer*exp(zs*res0_sdloer);
      
      List res1 = proc_OERvec(W_k, dd_i, dd_j);
      NumericVector res1_oer = res1["OER"];
      NumericVector res1_sdloer = res1["sdlOER"];
      NumericVector cc1_lb = res1_oer*exp(zs*res1_sdloer);
      
      for(int itt = 0 ; itt < (cc1_lb.size()-1); ++itt){
        int i = dd_i[itt];
        int j = dd_j[itt];
        if(lgW_0(i,j)>0){
          
          if(((lgW_1(i,j)==0)|((cc0_lb[itt] <1.0) & (cc1_lb[itt]<1.0)))) lgW_1(i,j) = 0;
          
          if( (cc0_lb[itt] <1.0) & (cc1_lb[itt]<1.0) ){
            rec_k.push_back(dd_k[1]+1);
            rec_i.push_back(i+1);
            rec_j.push_back(j+1);
          }
          
        }
      }
    }
    
    
  }
  
  DataFrame recs = DataFrame::create(_["k"]=rec_k,_["i"]=rec_i,_["j"]=rec_j);
  
  return List::create(_["lgW_0"]=lgW_0, _["lgW_1"]=lgW_1, _["CFrec"]=recs);
}

