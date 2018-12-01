#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;
// [[Rcpp::export]]
extern "C"{
  MatrixXd matMultC(MatrixXd A, MatrixXd B)
  {
    MatrixXd C=A*B;
    return C;
  }
  
}
