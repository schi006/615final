#include <Rcpp.h>
using namespace Rcpp;
void quickSort(NumericVector A, NumericVector B, int pl, int r);
NumericVector orderC(NumericVector A);
NumericMatrix cbindC(NumericMatrix mf, NumericMatrix Y);
//[[Rcpp::export]]
NumericMatrix cbindC(NumericMatrix mf, NumericMatrix Y) {
  int acoln = mf.ncol();
  NumericMatrix out = no_init_matrix(a.nrow(), acoln + 2);
  for (int j = 0; j < acoln + 2; j++) {
    if (j < 2) {
      out(_, j) = Y(_, 1-j);
    } else {
      out(_, j) = b(_, j - 2);
    }
  }
  return out;
}
NumericVector orderC(NumericVector A){
 int n=A.length();
 NumericVector B(n);
 for(int i=0;i<n;i++) B[i]=i+1;
 quickSort(A,B,0,n-1); 
 return B;
}
void quickSort(NumericVector A, NumericVector B, int pl, int r){
 if(pl<r){
  double piv=A[r];
  double pivB=B[r];
  int i=pl-1;
  double tmp,tmpB;
  for(int j=pl;j<r;++j){
   if(A[j]<piv){
    ++i;
    tmp=A[i];A[i]=A[j];A[j]=tmp;
    tmpB=B[i];B[i]=B[j];B[j]=tmpB;
   }
  }
  A[r]=A[i+1];A[i+1]=piv;
  B[r]=B[i+1];B[i+1]=pivB;
  quickSort(A,B,pl,i);
  quickSort(A,B,i+2,r);
 }
}
