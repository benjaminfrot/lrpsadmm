# lrpsadmm
Low-rank plus sparse estimation via Alternating Direction Method of Multipliers (ADMM)

An R package to fit the estimator suggested in Chandrasekaran et al. (https://projecteuclid.org/euclid.aos/1351602527). 
See the [manual](https://github.com/benjaminfrot/lrpsadmm/blob/master/manual.pdf) or the [lrpsadmm-examples](https://github.com/benjaminfrot/lrpsadmm-examples) for examples and more details. 

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

## Tutorials and Examples
See the [lrpsadmm-examples](https://github.com/benjaminfrot/lrpsadmm-examples) repository for step by step tutorials on how to use the package. The [manual](https://github.com/benjaminfrot/lrpsadmm/blob/master/manual.pdf) also contains examples.

## Description
Fit the Low-Rank plus Sparse (LRpS) estimator using the Alternating 
Direction Method of Multipliers (ADMM). This model learns an inverse 
covariance matrix which is the sum of a sparse matrix and a low-rank matrix as suggested by 
Chandrasekaran et al (2012) (DOI:10.1214/11-AOS949). The package supports robust estimation via an 
estimator of the correlation matrix based on Kendall rank correlations. It includes functions to compute 
a whole regularisation path and to select tuning parameters with cross-validation.

## Remarks and Issues
For questions and bug reporting, feel free to either send an email to benjamin dot frot @ gmail dot com or
open an issue.
