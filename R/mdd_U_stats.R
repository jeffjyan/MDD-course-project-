#' Unbiased martingale difference divergence statistics
#'
#' @description Unbiased estimate of squared martingale difference divergence and bias-corrected estimate of squared martingale difference correlation
#' @param X distances or data of first sample
#' @param Y distances or data of second sample
#'
#' @return \code{mdd_U} returns a list containing
#' \itemize{
#'   \item{MDD_sq_U: }{unbiased estimate of squared martingale difference divergence}
#'   \item{MDC_sq_U: }{bias-corrected estimate of squared martingale difference correlation}
#' }
#' @export
#'
#' @examples
#' A <- iris[1:50, 1:2]
#' B <- iris[1:50, 3:4]
#' mdd_U(A, B)
#' 
#' @references Park, Trevor, Xiaofeng Shao, and Shun Yao. "Partial martingale difference correlation." Electronic Journal of Statistics 9.1 (2015): 1492-1517.
mdd_U <- function(X, Y){
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
  
  # U-centered version of X and Y
  X_center <- U_center(X)
  Y <- Y^2 / 2
  Y_center <- U_center(Y)
  
  # Calculate unbiased martingale difference divergence statistics
  MDD_sq_U <- U_innerproduct(X_center, Y_center)
  X_norm <- sqrt(U_innerproduct(X_center))
  Y_norm <- sqrt(U_innerproduct(Y_center))
  if (X_norm * Y_norm == 0){
    MDC_sq_U <- 0
  }else {
    MDC_sq_U <- MDD_sq_U / (X_norm * Y_norm)
  }
  
  
  # Return unbiased martingale difference divergence statistics
  return(list(MDD_sq_U = MDD_sq_U, MDC_sq_U = MDC_sq_U))
}