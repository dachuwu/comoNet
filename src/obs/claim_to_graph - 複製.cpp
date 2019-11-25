#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List iseq_els(IntegerVector t, IntegerVector dz, bool drop_bilink){
  List els_i;
  if(drop_bilink){
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
  
  return(els_i);
}

// [[Rcpp::export]]
IntegerMatrix append_els_to_A(IntegerMatrix A0, List els){
  IntegerMatrix A1( A0.nrow() );
  for(int i=0; i< els.length(); ++i){
    
    List b = els[i];
    int r = b[0];
    IntegerVector cc = b[1];
    
    for(int j=0; j < cc.length(); ++j){
      A1((r-1),(cc[j]-1)) += 1 ;
    }
  }
  A1 += A0;
  return(A1);
}


// [[Rcpp::export]]
IntegerMatrix multiEls_to_AA(List els_ls, int Nc){
  IntegerMatrix A1 (Nc);
  for(int i=0; i< els_ls.length(); ++i){
    List els = els_ls[i];
    A1 = append_els_to_A(A1, els);
  }
  return(A1);
}
