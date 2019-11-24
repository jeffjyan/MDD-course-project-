#' Double centering
#' 
#' @description Double centered version of distance matrix of A
#' @param A dist object or data matrix
#'
#' @return \code{D_center} returns double centered distance matrix
#' @export
#'
#' @examples 
#' A <- iris[1:50, 1:4]
#' A_Dcenter <- D_center(A)
#' 
#' @references Park, Trevor, Xiaofeng Shao, and Shun Yao. "Partial martingale difference correlation." Electronic Journal of Statistics 9.1 (2015): 1492-1517.
D_center <- function(A){
  # Check whether A is a dist object or data matrix
  if (class(A) != "dist") A <- dist(A)
  A <- as.matrix(A)
  
  # grand mean, row mean and column mean
  A_mean <- mean(A)
  A_mean_r <- rowMeans(A)
  A_mean_c <- colMeans(A)
  
  # Return double centered distance matrix
  return(t(t(A) - A_mean_c) - A_mean_r + A_mean)
}


#' U-centering
#' 
#' @description U-centered version of distance matrix of A
#' @param A dist object or data matrix
#'
#' @return \code{U_center} returns U-centered distance matrix
#' @export
#'
#' @examples
#' A <- iris[1:50, 1:4]
#' A_Ucenter <- U_center(A)
#' 
#' @references Park, Trevor, Xiaofeng Shao, and Shun Yao. "Partial martingale difference correlation." Electronic Journal of Statistics 9.1 (2015): 1492-1517.
U_center <- function(A){
  # Check whether A is a dist object or data matrix
  if (class(A) != "dist") A <- dist(A)
  A <- as.matrix(A)
  
  n_data <- nrow(A)
  
  # grand sum, row sum and column sum
  A_sum <- sum(A)
  A_sum_r <- rowSums(A)
  A_sum_c <- colSums(A)
  
  # Calculate and return U-centering matrix
  A_U <- t(t(A) - A_sum_c/(n_data-2)) - A_sum_r/(n_data-2) + A_sum/((n_data-1)*(n_data-2))
  diag(A_U) <- 0
  return(A_U)
}