// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace std;
// [[Rcpp::export]]
arma::vec timesTwo(int x) {
  arma::vec t = arma::zeros(5);
  t += 2;
  t(3) = 3;
  cout << t(arma::span(0,-1)) <<endl;
  cout << max(t) <<endl;
  if( true & true ){
    cout << "max(t)" <<endl;
  }
  for( int k=0; k< 0; ++k){
    cout << k << endl;
  }
  
  //return Rcpp::runif(5,0,0.5);
  return Rcpp::runif(2,0,0.5);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
