mdd <- function(X, Y){
  # compatibility check
  if(nrow(X) != nrow(Y)){
    stop("Sample sizes must agree")
  }
  if (!(all(is.finite(X))) | !(all(is.finite(Y)))){
    stop("Data contains missing or infinite values")
  }
  
  n <- nrow(X)
  
  Dcenter <- function(X){
    # grand mean, row mean and column mean
    X_mean <- mean(X)
    X_mean_r <- rowMeans(X)
    X_mean_c <- colMeans(X)
    
    # return double centering matrix
    return(t(t(X) - X_mean_c) - X_mean_r + X_mean)
  }
  
  Xd <- dist(Xd)
  Xd_center <- Dcenter(Xd)
  Yd <- dist(Yd)
  Yd <- Yd^2 / 2
  Yd_center <- Dcenter(Y)
  
  MDD <- sqrt(mean(Xd * Yd))
  dVarX <- sqrt(mean(Xd^2))
  dVarY <- sqrt(mean(Yd^2))
  dVarXY <- sqrt(dVarX * dVarY)
  if(dVarXY > 0){
    MDC <- MDD / dVarXY
  }else{
    MDC <- 0
  }
  
  return(list(MDD = MDD, MDC = MDC))
}