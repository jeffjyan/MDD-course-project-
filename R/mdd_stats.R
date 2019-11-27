#' Martingale Difference Divergence and Correlation Statistics
#' 
#' @description Sample martingale difference divergence and sample martingale difference correlation between X and Y
#' @param X distances or data of first sample
#' @param Y distances or data of second sample
#'
#' @return \code{mdd} returns a list containing
#' \itemize{
#'   \item{MDD: }{sample martingale difference divergence}
#'   \item{MDC: }{sample martingale difference correlation}
#' }
#' @export
#'
#' @examples
#' A <- iris[1:50, 1:2]
#' B <- iris[1:50, 3:4]
#' mdd(A, B)
#' 
#' @references Shao, Xiaofeng, and Jingsi Zhang. "Martingale difference correlation and its use in high-dimensional variable screening." Journal of the American Statistical Association 109.507 (2014): 1302-1318.
mdd <- function(X, Y){
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
  
  # double centered version of X and Y
  X_center <- D_center(X)
  Y <- Y^2 / 2
  Y_center <- D_center(Y)
  
  # Calculate the sample martingale difference divergence and correlation
  MDD <- sqrt(mean(X_center * Y_center))
  dVarX <- sqrt(mean(X_center^2))
  VarY <- sqrt(mean(Y_center^2))
  dVarXY <- sqrt(dVarX * VarY)
  if (dVarXY > 0){
    MDC <- MDD / dVarXY
  }else if (dVarXY == 0){
    MDC <- 0
  }
  
  # Return sample martingale difference divergence and correlation
  return(list(MDD = MDD, MDC = MDC))
}