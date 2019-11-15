Dcenter <- function(X){
  # grand mean, row mean and column mean
  X_mean <- mean(X)
  X_mean_r <- rowMeans(X)
  X_mean_c <- colMeans(X)
  
  # return double centering matrix
  return(t(t(X) - X_mean_c) - X_mean_r + X_mean)
}

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