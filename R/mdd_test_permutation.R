#' Conditional mean independence test using permutation
#'
#' @param X distances or data of first sample
#' @param Y distances or data of second sample
#' @param R number of replicates
#'
#' @return \code{mdd_test_perm} returns a list containing
#' \itemize{
#'   \item{test_statistic: }{observed value of the test statistic}
#'   \item{MDD: }{sample martingale difference divergence}
#'   \item{replicates: }{replicates of the test statistic}
#'   \item{p_value: }{approximate p-value of the test}
#' }
#' @export
#'
#' @examples
#' A <- iris[1:50, 1:2]
#' B <- iris[1:50, 3:4]
#' set.seed(1)
#' mdd_test_perm(A, B)
#' 
#' @references Shao, Xiaofeng, and Jingsi Zhang. "Martingale difference correlation and its use in high-dimensional variable screening." Journal of the American Statistical Association 109.507 (2014): 1302-1318.
#' @references Park, Trevor, Xiaofeng Shao, and Shun Yao. "Partial martingale difference correlation." Electronic Journal of Statistics 9.1 (2015): 1492-1517.
mdd_test_perm <- function(X, Y, R = 999){
  # Check whether X and Y are dist objects or data matrices
  if (!(class(X) == "dist")) X <- dist(X)
  if (!(class(Y) == "dist")) Y <- dist(Y)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # compatibility check
  if (nrow(X) != nrow(Y)){
    stop("Sample sizes must agree")
  }
  if (!(all(is.finite(X))) | !(all(is.finite(Y)))){
    stop("Data contains missing or infinite values")
  }
  
  # double centering
  X_center <- D_center(X)
  Y <- Y^2 / 2
  Y_center <- D_center(Y)
  
  # Conduct the test using approxiamte permutation test
  return(mdd_test_perm_c(X_center, Y_center, R))
}