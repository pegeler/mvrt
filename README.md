# mvrt: An R package to generate multivariate t-distributed random data.

You can see an exploration of how this package is used by reading the [package vignette](https://pegeler.github.io/mvrt/). Different methods for function definitions are used and performance is assessed.

## Example Usage

```r
S <- convert_R2S(make_cor_mat(.9),2:3)
x <- mvrt( 30, 4:5, S, 29)
y <- mvrt2(30, 4:5, S, 29, .01)

cor(x)
##           [,1]      [,2]
## [1,] 1.0000000 0.8884765
## [2,] 0.8884765 1.0000000

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

