#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// conditional mean independence test based on martingale difference divergence statistics, using permutation test
//
// [[Rcpp::export]]
Rcpp::List mdd_test_perm_c(const arma::mat& X, const arma::mat& Y, int R) {
  double p_value;
  arma::vec replicates(R);
  
  // number of data
  int n = X.n_rows;
  
  // Calculate sample martingale difference divergence
  double MDD = sqrt(accu(X % Y)) / n;
  
  // Calculate test statistic
  double test_stat = accu(X % Y) / n;
  
  // Perform approxiamte permutation test
  for(int r = 0; r < R; r++){
    double sum = 0;
    
    // Generate a random permutation
    arma::uvec gen_ind = arma::randperm(n);
    
    // Calculate test statistic for each permutation
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++){
        sum += X(gen_ind[i],gen_ind[j]) * Y(i,j);
      }
    replicates[r] = sum / n;
  }
  
  // Calcualte p-value
  p_value = (1.0 + accu(replicates >= test_stat)) / (1.0 + R);
  
  return List::create(Rcpp::Named("test_statistic") = test_stat,
                      Rcpp::Named("MDD") = MDD, 
                      Rcpp::Named("replicates") = replicates, 
                      Rcpp::Named("p_value") = p_value);
}