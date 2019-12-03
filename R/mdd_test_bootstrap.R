#' Conditional mean independence test using wild bootstrap
#'
#' @param X distances or data of first sample
#' @param Y distances or data of second sample
#' @param B number of replicates
#'
#' @return \code{mdd_test_boot} returns a list containing
#' \itemize{
#'   \item{test_statistic: }{observed value of the test statistic}
#'   \item{MDD_sq_U: }{unbiased estimate of squared martingale difference divergence}
#'   \item{replicates: }{replicates of the test statistic}
#'   \item{p_value: }{approximate p-value of the test}
#' }
#' @export
#'
#' @examples
#' A <- iris[1:50, 1:2]
#' B <- iris[1:50, 3:4]
#' set.seed(1)
#' mdd_test_boot(A, B)
#' 
#' @references Chung Eun Lee, Xianyang Zhang, and Xiaofeng Shao. Testing the conditional mean independence for functional data. Biometrika, forthcoming, 2019.
mdd_test_boot <- function(X, Y, B = 999){
  # Check whether X and Y are dist objects or data matrices
  if (!(class(X) == "dist")) X <- dist(X)
  if (!(class(Y) == "dist")) Y <- dist(Y)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # compatibility check
  if(nrow(X) != nrow(Y)){
    stop("Sample sizes must agree")
  }
  if (!(all(is.finite(X))) | !(all(is.finite(Y)))){
    stop("Data contains missing or infinite values")
  }
  
  # U-centering
  X_center <- U_center(X)
  Y <- Y^2 / 2
  Y_center <- U_center(Y)
  
  # Conduct the test using wild Bootstrap
  return(mdd_test_boot_c(X_center, Y_center, B))
}