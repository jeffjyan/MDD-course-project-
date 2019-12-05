# MDD: Conditional Mean Independence Testing via Martingale Difference Divergence

The R package `MDD` implements conditional mean independence testing based on martingale difference divergence, using permutation and wild Bootstrap. 

## Introduction and Motivation
Martingale difference divergence (MDD) and martingale dfference correlation (MDC) are intended to measure departure from the relationship *E(Y|X) = E(Y) a.s.*. For example, in linear models, an important problem is to assess whether *X* contributes to the conditional mean of *Y*, i.e., *H0 : E(Y|X) = E(Y)*. Since if *H0* is supported by the data, then there is no need to pursue a regression model for the mean of *Y* given *X*.

## Installation
```r
require(devtools)
devtools::install_github("jeffjyan/MDD")
library(MDD)
```

## Example
```r
# Data illustration
X <- mtcars$wt
Y <- mtcars$mpg
mdd_test_perm(X, Y)$p_value
mdd_test_boot(X, Y)$p_value

X <- mtcars$am
Y <- mtcars$carb
mdd_test_perm(X, Y)$p_value
mdd_test_boot(X, Y)$p_value

X <- iris[1:50, c(1,3)]
Y <- iris[1:50, c(2,4)]
mdd_test_perm(X, Y)$p_value
mdd_test_boot(X, Y)$p_value


# Simulation
X <- rnorm(30)
Y <- rpois(30, lambda = 1)
mdd_test_perm(X, Y)$p_value
mdd_test_boot(X, Y)$p_value

X <- rnorm(30)
Y <- X^2
mdd_test_perm(X, Y)$p_value
mdd_test_boot(X, Y)$p_value
```

## Reference
Székely, Gábor J., Maria L. Rizzo, and Nail K. Bakirov. "Measuring and testing dependence by correlation of distances." The annals of statistics 35.6 (2007): 2769-2794. 

Shao, Xiaofeng, and Jingsi Zhang. "Martingale difference correlation and its use in high-dimensional variable screening." Journal of the American Statistical Association 109.507 (2014): 1302-1318.

Park, Trevor, Xiaofeng Shao, and Shun Yao. "Partial martingale difference correlation." Electronic Journal of Statistics 9.1 (2015): 1492-1517.

Chung Eun Lee, Xianyang Zhang, and Xiaofeng Shao. Testing the conditional mean independence for functional data. Biometrika, forthcoming, 2019.