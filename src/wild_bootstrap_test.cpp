#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// conditional mean independence test based on unbiased martingale difference divergence statistics, using wild bootstrap
//
// [[Rcpp::export]]
Rcpp::List mdd_test_boot_c(const arma::mat& X, const arma::mat& Y, int B) {
  double p_value;
  arma::vec replicates(B);
  
  // number of data
  int n = X.n_rows;
  
  // Calculate unbiased estimate of squared martingale difference divergence
  double MDD_sq_U = accu(X % Y) / (n * (n-3));
  
  // Calculate test statistic
  double test_stat = MDD_sq_U * n;
  
  // Perform wild Bootstrap test
  for(int b = 0; b < B; b++){
    double sum = 0;
    
    // Generate a sequence of standard normal distributions
    arma::vec gen_norm = arma::randu(n);
    
    // Calculate test statistic
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++){
        sum += gen_norm[i] * X(i,j) * Y(i,j) * gen_norm[j];
      }
    replicates[b] = sum / (n - 3);
  }
  
  // Calcualte p-value
  p_value = (1.0 + accu(replicates >= test_stat)) / (1.0 + B);
  
  return List::create(Rcpp::Named("test_statistic") = test_stat,
                      Rcpp::Named("MDD_sq_U") = MDD_sq_U, 
                      Rcpp::Named("replicates") = replicates, 
                      Rcpp::Named("p_value") = p_value);
}