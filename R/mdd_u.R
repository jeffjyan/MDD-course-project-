mdd_u <- function(X, Y){
  # compatibility check
  if(nrow(X) != nrow(Y)){
    stop("Sample sizes must agree")
  }
  if (!(all(is.finite(X))) | !(all(is.finite(Y)))){
    stop("Data contains missing or infinite values")
  }
  
  n <- nrow(X)
  
  Ucenter <- function(X){
    X_n <- nrow(X)
    
    # grand sum, row sum and column sum
    X_sum <- sum(X)
    X_sum_r <- rowSums(X)
    X_sum_c <- colSums(X)
    
    # calculate U-centering matrix
    Xd <- t(t(X) - X_sum_c/(X_n-2)) - X_sum_r/(X_n-2) + X_sum/((X_n-1)*(X_n-2))
    diag(Xd) <- 0
    return(Xd)
  }
  
  Xd <- dist(Xd)
  Xd_center <- Ucenter(Xd)
  Yd <- dist(Yd)
  Yd <- Yd^2 / 2
  Yd_center <- Ucenter(Y)
  
  MDD_u <- sqrt(sum(Xd * Yd) / (n*(n-3)))
  
  return(MDD_u)
}