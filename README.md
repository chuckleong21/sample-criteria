<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/chuckleong21/sample-criteria/graph/badge.svg)](https://app.codecov.io/gh/chuckleong21/sample-criteria)
  <!-- badges: end -->

# sample.criteria
The goal of `sample.criteria` is to calculate the minimum sample size based on reported clinic prediction models for future work.

## Installaion
The package is not on CRAN, install it through github:

```r
remotes::install_github("jingwen517/sample-criteria")
```

## Usage 
To calculate the minimum sample provided only C statistics and outcome proportion are reported for a binary outcome logistic regression model: 

```r
library(sample.criteria)

set.seed(1234)
pmsamplesize(Q = 30, k = 2, auc = 0.81, prev = 0.77)
```

To learn more about other usage cases, see `vignette("sample-criteria")`.
