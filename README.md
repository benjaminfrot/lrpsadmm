# lrpsadmm
Low-rank plus sparse estimation via Alternating Direction Method of Multipliers (ADMM)

An R package to fit the estimator suggested in Chandrasekaran et al. (https://projecteuclid.org/euclid.aos/1351602527). 
See the [manual](https://github.com/benjaminfrot/lrpsadmm/blob/master/manual.pdf) for examples and more details.

## Installation

(Optional) Install the `devtools` package:
```R
install.packages("devtools")
```
Install the package:
```R
library(devtools)
install_github("benjaminfrot/lrpsadmm")
```

## Tutorials
See this repository : https://github.com/benjaminfrot/lrpsadmm-examples

## Description
Fit the Low-Rank plus Sparse (LRpS) estimator using an accelerated version of the Alternating 
Direction Method of Multipliers (ADMM) (DOI:10.1137/120896219). This model learns an inverse 
covariance matrix which is the sum of a sparse matrix and a low-rank matrix as suggested by 
Chandrasekaran et al (2012) (DOI:10.1214/11-AOS949). The package supports robust estimation via an 
estimator of the correlation matrix based on Kendall rank correlations. It includes functions to compute 
a whole regularisation path and to select tuning parameters with cross-validation.

## Remarks and Issues
For questions and bug reporting, feel free to either send an email to benjamin dot frot @ gmail dot com or
open an issue.
