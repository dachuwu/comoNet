#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List filter_OER(arma::mat A, double sig_lv){
  mat OER= mat(A.n_rows, A.n_cols, fill::zeros);
  mat varlOER= mat(A.n_rows, A.n_cols, fill::zeros);
  int nn = A.n_rows;
  double S = as_scalar(sum(sum(A, 1)));
  colvec s_i = sum(A, 1);
  rowvec s_j = sum(A, 0);
  for(int i=0; i < nn; ++i){
    for(int j=0; j < nn; ++j){
      double a = A(i,j)/S;
      if(a==0) continue;
      double b = s_i(i)/S - a;
      double c = s_j(j)/S - a;
      OER(i,j) = a/((a+b)*(a+c));
      varlOER(i,j) = (b*c*(1-a) + (1-b-c)*pow(a,2) - pow(a,3))/(a*(a+b)*(a+c)*S);
    }
  }
  double z =  R::qnorm(sig_lv/2, 0.0, 1.0, 1, 0);
  mat lb =  OER % exp(z*sqrt(varlOER)) ;
  umat logiB =  lb > ones(A.n_rows, A.n_cols);
  mat B =  A % conv_to<dmat>::from(logiB);
  
  return List::create(_["Raw_OER"]=OER, _["A_filtered"]=B);
}

double Ftest_c(arma::mat x) {
  Function ft("fisher.test");   
  List res = ft(x,_["alternative"]="g",_["conf.int"]=false);
  return res["p.value"];
}

// [[Rcpp::export]]
List filter_phi(arma::mat A, double sig_lv){
  mat Fpval= ones(A.n_rows, A.n_cols) ;
  mat Raw_phi= zeros(A.n_rows, A.n_cols) ;
  int nn = A.n_rows;
  double S = as_scalar(sum(sum(A, 1)));
  colvec s_i = sum(A, 1);
  rowvec s_j = sum(A, 0);
  for(int i=0; i < nn; ++i){
    for(int j=0; j < nn; ++j){
      double a = A(i,j);
      if(a==0) continue;
      mat x ;
      x << a << s_i(i)-a << endr
        << s_j(j)-a << S-s_i(i)-s_j(j)+a << endr;
      Fpval(i,j) = Ftest_c(x);
      Raw_phi(i,j) = (a*S + s_i(i)*s_j(j))/sqrt(s_i(i)*s_j(j)*(S-s_i(i))*(S-s_j(j)));
    }
  }
  umat logiB =  Fpval < sig_lv;
  mat B =  A % conv_to<dmat>::from(logiB);
  return List::create(_["Raw_phi"] = Raw_phi, 
                      _["A_filtered"] = B,
                      _["Ft_pval"] = Fpval);
}

// [[Rcpp::export]]
List filter_dispar(arma::mat A, double sig_lv){
  mat Dspval= ones(A.n_rows, A.n_cols) ;
  
  int nn = A.n_rows;
  umat logiA = A>0;
  colvec s_i = sum(A, 1);
  colvec k_i = sum(conv_to<dmat>::from(logiA), 1);
  for(int i=0; i < nn; ++i){
    for(int j=0; j < nn; ++j){
      double a = A(i,j);
      if(a==0) continue;
      Dspval(i,j) = pow(1-a/s_i(i), k_i(i)-1);
    }
  }
  umat logiB =  Dspval < sig_lv;
  mat B =  A % conv_to<dmat>::from(logiB);
  return List::create(_["A_filtered"] = B,
                      _["Ds_pval"] = Dspval);
}
