## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(sample.criteria)

## -----------------------------------------------------------------------------
# remotes::install_github("your.name/your.repository")

## -----------------------------------------------------------------------------
library(sample.criteria)
# Noted here that k = 2 for binary outcomes
approximate_R2(k = 2, auc = 0.81, prev = 0.77)

## -----------------------------------------------------------------------------
set.seed(1234)
scriteria <- pmsamplesize(Q = 30, k = 2, auc = 0.81, prev = 0.77, r2_nagelkerke = 0.15)
scriteria 

## -----------------------------------------------------------------------------
pmsamplesize(Q = 24, p = c(138 - 24, 24), r2_nagelkerke = 0.48)

## -----------------------------------------------------------------------------
# Note the prev argument is not default to NULL
set.seed(101)
a <- pmsamplesize(Q = 17,
             k = 5,
             p = c(0.729, 0.053, 0.05, 0.133, 0.034),
             r2_nagelkerke = 0.15,
             shrinkage = 0.9,
             auc = c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82))

