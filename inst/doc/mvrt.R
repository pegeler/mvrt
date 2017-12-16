## ---- include=FALSE------------------------------------------------------
knitr::opts_chunk$set(fig.width = 6, fig.height = 4, fig.align = 'center')

## ------------------------------------------------------------------------
n <- 30
mu <- c(5,4)
R <- matrix(c(1,.9,.9,1),2,2)
var <- c(2.5,2)

S <- mvrt::convert_R2S(R, var)

## ------------------------------------------------------------------------
# With a for loop
mvrt_R_a <- function(n, mu, S, df = n - 1) {

  g.t <- t(chol(S))

  bivMat <- matrix(0,nrow = length(mu), ncol = n)

  for (i in seq_len(n)) {
    bivMat[,i] <- mu + g.t %*% rt(length(mu), df)
  }

  t(bivMat)

}

# Without a for loop
mvrt_R_b <- function(n, mu, S, df = n - 1) {
  
  g <- chol(S)
  
  random_matrix <- matrix(rt(n*length(mu), df),nrow = length(mu))
  deviation <- t(g) %*% random_matrix
  
  t(mu + deviation)
  
}

## ---- results='hide', message=FALSE--------------------------------------
Rcpp::sourceCpp("../src/mvrt.cpp")

## ---- echo=FALSE, results="asis"-----------------------------------------
cat("```c\n")
cat(readLines("../src/mvrt.cpp"), sep = "\n")
cat("```\n")

## ------------------------------------------------------------------------
get_cor_dist <- function(FUN, times, ..., seed = 999) {
  
  set.seed(seed)
  
  FUN <- match.fun(FUN)
  
  my_call <- as.call(list(FUN, ...))
  
  out <- numeric(times)
  
  for (i in seq_len(times)) {
    
    out[i] <- cor(eval(my_call))[2]
    
  }
  
  out
  
}

## ------------------------------------------------------------------------
mvrtRa_dist <-  get_cor_dist(mvrt_R_a, 100, n, mu, S)
mvrtRb_dist <-  get_cor_dist(mvrt_R_b, 100, n, mu, S)
mvrt_dist <-    get_cor_dist(mvrt, 100, n, mu, S, n - 1)  
MASS_dist <-    get_cor_dist(MASS::mvrnorm, 100, n, mu, S)
mvtnorm_dist <- get_cor_dist(mvtnorm::rmvt, 100, n, S, n - 1)

## ------------------------------------------------------------------------
if (identical(mvrtRa_dist, mvrt_dist) && identical(mvrtRb_dist, mvrt_dist)) {
  message("Distributions of correlation coefficients are identical")
} else {
  message("Distributions of correlation coefficients are NOT identical")
}

if (
  identical(
    (mvrt_sample <- {set.seed(999); mvrt(n, mu, S, n - 1)}),
    {set.seed(999); mvrt_R_a(n, mu, S)}
    ) &&
  identical(
    mvrt_sample,
    {set.seed(999); mvrt_R_b(n, mu, S)}
    )
) {
  message("Random samples generated match")
} else {
  message("Random samples generated do NOT match")
}

## ------------------------------------------------------------------------
# Mean mu
colMeans(t(sapply(lapply(rep(30,100), mvrt, mu, S, n - 1),colMeans)))

# Mean var
colMeans(t(sapply(lapply(rep(30,100), mvrt, mu, S, n - 1), apply, 2, var)))

## ------------------------------------------------------------------------
mvrt_sample <- mvrt(n, mu, S, n - 1)

colMeans(mvrt_sample)
apply(mvrt_sample, 2,var)
cor(mvrt_sample)

## ------------------------------------------------------------------------
plot(mvrt_sample)
# Regression line
abline(lm(mvrt_sample[,2] ~ mvrt_sample[,1]), col = "red")

# Center of points
points(mean(mvrt_sample[,1]), mean(mvrt_sample[,2]), col = "blue", pch = 22)

## ------------------------------------------------------------------------
# Summary of the distributions
summary(mvrt_dist)
summary(MASS_dist)
summary(mvtnorm_dist)

# Graphical representations
my_hist <- function (data) {
  hist(
    data, 
    xlim = c(-1,1), 
    breaks = "fd", 
    main = paste("Histogram of",deparse(substitute(data))),
    xlab = deparse(substitute(data))
  )
  
  abline(v = 0.9, col = 'red')
}

my_hist(mvrt_dist)
my_hist(MASS_dist)
my_hist(mvtnorm_dist)

## ------------------------------------------------------------------------
# Define a function for making a 2x2 correlation matrix
make_cor_mat <- function (r) {  
  out <- matrix(r, 2, 2)
  diag(out) <- 1
  out
}

# A single reference matrix
R_ref <- make_cor_mat(.9)

# A set of test matrices
R_test <- lapply(seq(.5,1,.05), make_cor_mat)

# What does the determinant look like for the difference matrices?
sapply(R_test, function(a,b) det(a - b), R_ref)

# And the matrix norm?
sapply(R_test, function(a,b) norm(a - b, "m"), R_ref)

## ------------------------------------------------------------------------
# Using recursive function calls
mvrt_R_c <- function(n, mu, S, df = n - 1, max.norm = 0.05, type.norm = "m") {
 
  g <- chol(S)
 
  random_matrix <- matrix(rt(n*length(mu), df),nrow = length(mu))
  deviation <- t(g) %*% random_matrix
 
  out <- t(mu + deviation)
 
  if (norm(cov2cor(S) - cor(out), type = type.norm) <= max.norm) {
    return(out)
  } else {
    eval(match.call())
  }
 
}


# Using the while loop
mvrt_R_d <- function(n, mu, S, df = n - 1, max.norm = 0.05, type.norm = "m") {
  
  g.t <- t(chol(S))
  
  accept_data <- FALSE
  
  R_ref <- cov2cor(S)
  
  while (!accept_data) {
    
    random_matrix <- matrix(rt(n*length(mu), df),nrow = length(mu))
    deviation <- g.t %*% random_matrix
    
    out <- t(mu + deviation)
    
    accept_data <- norm(R_ref - cor(out), type = type.norm) <= max.norm
      
  } 
  
  out
  
}


## ------------------------------------------------------------------------
mvrtRc_dist <-  get_cor_dist(mvrt_R_c, 100, n, mu, S)
mvrtRd_dist <-  get_cor_dist(mvrt_R_d, 100, n, mu, S)

if (identical(mvrtRc_dist, mvrtRd_dist)) {
  message("Distributions of correlation coefficients are identical")
} else {
  message("Distributions of correlation coefficients are NOT identical")
}

summary(mvrtRc_dist)

my_hist(mvrtRc_dist)

## ------------------------------------------------------------------------
summary(get_cor_dist(mvrt_R_d, 100, n, mu, mvrt::convert_R2S(make_cor_mat(0.2),var)))

## ------------------------------------------------------------------------
microbenchmark::microbenchmark(
  mvrt_R_a(n, mu, S),               # 'for' loop
  mvrt_R_b(n, mu, S),               # Optimized R code
  mvrt_R_c(n, mu, S),               # Recursive call
  mvrt_R_d(n, mu, S),               # 'while' loop
  mvrt(n, mu, S, n - 1),            # C++ code
  MASS::mvrnorm(n, mu, S),
  mvtnorm::rmvt(n, S, n - 1) + mu
)

## ------------------------------------------------------------------------
MASS_mvrnorm <- MASS::mvrnorm
mvtnorm_rmvt <- mvtnorm::rmvt

microbenchmark::microbenchmark(
  MASS::mvrnorm(n, mu, S),
  MASS_mvrnorm(n, mu, S),
  mvtnorm::rmvt(n, S, n - 1) + mu,
  mvtnorm_rmvt(n, S, n - 1) + mu
)

