#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

List iseq_els(IntegerVector t, IntegerVector dz, bool bidir){
  List els_i;
  if(!bidir){
    for(int i=0; i < t.length()-1; ++i){
      int X = dz[i];
      IntegerVector Y = dz[t > t[i]];
      Y = Y[Y != dz[i]];
      List out = List::create(X, Y);
      els_i.push_back(out);
    }
  } else {
    for(int i=0; i < t.length()-1; ++i){
      int X = dz[i];
      IntegerVector Y = dz[t >= t[i]];
      Y = Y[Y != dz[i]];
      List out = List::create(X, Y);
      els_i.push_back(out);
    }
  }
  
  return els_i;
}


umat append_els_to_A(umat A0, List els){
  umat A1(A0.n_rows, A0.n_cols, fill::zeros);
  
  for(int l=0; l< els.length(); ++l){
    
    List b = els[l];
    int r = b[0];
    IntegerVector cc = b[1];
    
    for(int j=0; j < cc.length(); ++j){
      A1((r-1),(cc[j]-1)) += 1 ;
    }
  }
  A1 = A1 + A0;
  return A1;
}

// [[Rcpp::export]]
List build_Raw(List iseq, StringVector dzlv, bool bidir=false){

  int nn = dzlv.length();
  umat W(nn, nn, fill::zeros);
  
  for(int i=0; i< iseq.length(); ++i){
    DataFrame tdz_i = iseq[i];
    if(tdz_i.nrows() <=1) continue;
    W= append_els_to_A(W, iseq_els(tdz_i["t"], tdz_i["dz"], bidir));
  }
  return  List::create(W);
}

// [[Rcpp::export]]
List build_CondRaw(List iseq, StringVector dzlv, bool bidir=false){
  
  int nn = dzlv.length();
  umat W(nn, nn, fill::zeros);
  ucube U(nn, nn, nn, fill::zeros);
  
  for(int i=0; i< iseq.length(); ++i){
    
    DataFrame tdz_i = iseq[i];
    if(tdz_i.nrows() <=1) continue;
    W= append_els_to_A(W, iseq_els(tdz_i["t"], tdz_i["dz"], bidir));
    
    if(tdz_i.nrows() < 3) continue;
    for(int k=0; k< (tdz_i.nrows()-2); ++k){
      IntegerVector s_t = tdz_i["t"];
      IntegerVector s_dz =  tdz_i["dz"];
      int kdz = s_dz[k]-1;
      s_t.erase(0,k+1);
      s_dz.erase(0,k+1);
      U.slice(kdz) = append_els_to_A(U.slice(kdz), iseq_els(s_t, s_dz, bidir));
    }
  }
  
  return  List::create(W, U);
}
