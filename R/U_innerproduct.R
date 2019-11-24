#' Inner product in the Hilbert space of U-centered distance matrices
#'
#' @param A U-centered distance matrix
#' @param B U-centered distance matrix
#'
#' @return \code{U_innerproduct} returns the inner product of A and B, which is a scalar
#' @export
#'
#' @examples
#' A <- iris[1:50, 1:2]
#' B <- iris[1:50, 3:4]
#' A_Ucenter <- U_center(A)
#' B_Ucenter <- U_center(B)
#' U_innerproduct(A_Ucenter, B_Ucenter)
#' 
#' @references Park, Trevor, Xiaofeng Shao, and Shun Yao. "Partial martingale difference correlation." Electronic Journal of Statistics 9.1 (2015): 1492-1517.
U_innerproduct <- function(A, B = NULL){
  n_data <- nrow(A)
  
  # Check whether B is null or not. If B is null, then return the inner product of A and itself. 
  if (is.null(B)){
    return(sum(A^2) / (n_data * (n_data - 3)))
  }else{
    return(sum(A * B) / (n_data * (n_data - 3)))
  }
}