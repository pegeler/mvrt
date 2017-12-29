# mvrt: An R package to generate multivariate t-distributed random data.

## Introduction

This package provides random matrices with a defined covariance structure. The package began as an exercise to show performance improvements of using vectorized R code rather than explicit `for` loops. Later, it became a vehicle for working with the [`Rcpp`](https://cran.r-project.org/package=Rcpp) and [`RcppArmadillo`](https://cran.r-project.org/package=RcppArmadillo) packages and learning how to incorporate them into packages of my own. Future work may include `FORTRAN` code to further explore the capabilities of using compiled code in R scripts and packages.

You can see an exploration the thought process just described in the [package vignette](https://pegeler.github.io/mvrt/). Different methods for function definitions are used and performance is assessed.

## Example Usage

```r
library(mvrt)

S <- convert_R2S(make_cor_mat(.9),2:3) # Create a covariance matrix
x <- mvrt( 30, 4:5, S, 29)             # Generate random data with n = 30
y <- mvrt2(30, 4:5, S, 29, .01)        # Random data with maximum abs deviation
                                       # from input parameters specified

# Correlation matrix of x
cor(x)
##           [,1]      [,2]
## [1,] 1.0000000 0.8884765
## [2,] 0.8884765 1.0000000

# Correlation matrix of y
cor(y)
##          [,1]     [,2]
## [1,] 1.000000 0.905187
## [2,] 0.905187 1.000000
```

## Installation

I recommend that you have [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) installed in order to download and install this package. To do so, type the following into your R console:

```r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("pegeler/mvrt")
```

