#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
NumericMatrix cbindC(NumericMatrix mf, NumericMatrix Y);//Y has exactly 2 columns;
NumericVector matMultC(NumericMatrix aa, NumericVector bb);
NumericVector elementE(NumericVector x);
NumericMatrix nominatorC(NumericMatrix d, NumericVector s);
//[[Rcpp::export]]
NumericMatrix cbindC(NumericMatrix mf, NumericMatrix Y) {
  int acoln = mf.ncol();
  NumericMatrix out = no_init_matrix(Y.nrow(), acoln + 2);
  for (int j = 0; j < acoln + 2; j++) {
    if (j < 2) {
      out(_, j) = Y(_, 1-j);
    } else {
      out(_, j) = mf(_, j - 2);
    }
  }
  return out;
}
//[[Rcpp::export]]
NumericVector matMultC(NumericMatrix aa, NumericVector bb){
 Map<MatrixXd> am(as<Map<MatrixXd> >(aa));
 Map<VectorXd> bv(as<Map<VectorXd> >(bb));
 return(wrap(am*bv));
}
//[[Rcpp::export]]
NumericVector elementE(NumericVector x){
   int nrow=x.size();
   for(int i=0;i<nrow;i++){
         x(i)=exp(x(i));
   }
   return(x);
}
//[[Rcpp::export]]
NumericMatrix nominatorC(NumericMatrix d, NumericVector s) {
int coln = d.ncol();
int rown = d.nrow();
NumericMatrix out = no_init_matrix(rown, coln);
for(int k=0;k<coln;k++){
 out(0,k)=d(0,k)*s(0);
}
for(int i=1;i<rown;i++){
 for(int j=0;j<coln;j++){
  out(i,j)=d(i,j)*s(i)+out(i-1,j);
 }
}
return out;
}