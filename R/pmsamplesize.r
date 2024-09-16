#' @exportS3Method base::print
print.pmsample <- function(x, ...) { # nocov start
  message(paste("NB: Assuming", x$sigma, "acceptable difference in apparent & adjusted R-squared"))
  message(paste("NB: Assuming", x$margin_error, "margin of error in estimation of intercept"))
  message(paste("NB: Events per Predictor Parameter assumes prevalance", round(x$prev, 3), "\n"))
  print(x$criteria)
  cat(sprintf("Minimum sample size required for new model development based on user inputs = %.0f,
 with %.0f events (assuming an outcome prevalence = %.3g) and an EPP = %.3g\n",
              x$final$sample_size,
              x$final$sample_size * x$prev,
              x$prev, x$final$EPP))
} # nocov end

#' @exportS3Method base::print
print.pmsample_multi <- function(x, ...) { # nocov start
  message(paste("NB: Assuming", x$sigma, "acceptable difference in apparent & adjusted R-squared"))
  message(paste("NB: Assuming", x$margin_error, "margin of error in estimation of intercept"))
  print(x$criteria)
} # nocov end


#' Sample Criteria Extension for Multinomial Regression
#'
#' @importFrom Rdpack reprompt
#'
#' @description
#' Calculate sample criteria extended for multinominal regression.
#'
#' @details
#' The pmsample objects created by \code{pmsamplesize} can represent the minimum
#' sample size needed for developing a new model based on the current studies.
#'
#' \code{pmsamplesize} function adapts the amount of information based on
#' the type of model. A certain parameters need pre-specifying for calculation:
#' - The number of candidate predicators \emph{Q}, usually obtained using
#' the fractional polynomial method and it is not addressed by this package.
#' - The number of levels in response variable \emph{k}, typically for
#' a binary model k = 2, for multinomial models k > 2. One of the argument
#' between \code{k} or \code{p}, the outcome proportions or the observation counts
#' of the events, must be specified. In the case of a missing \code{k}, it will
#' then defined by the length of \code{p}.
#' - The order of p matters for calculation: for a binary model, the second value
#' in the vector is perceived as 'outcome proportion' or 'outcome events'.
#' - Pairwise C-statistics \code{auc} for simulation. The process requires
#' estimates of the pairwise outcome proportions \eqn{\phi_{k,r}} of
#' category \emph{k} relative to the reference category \emph{r}, \code{prev}.
#' If not supplied, they will be derived from \code{p}.
#'
#' @references
#' \insertCite{Pate2023-yh}{sample.criteria}
#' \insertCite{Riley2019-rn}{sample.criteria}
#'
#' @param Q An integer. The number of candidate predictors. See Details
#' @param k An integer. The number of levels in response variable
#' @param p A numeric vector. Outcome proportion(s) or event observations
#' @param adjust A logical. Should adjusted R-square value be used? Default to FALSE
#' @param r2_cs_app A numeric vector. The apparent Cox-Snell \eqn{R^2}
#' @param r2_cs_adj A numeric vector. The optimism-adjusted Cox-Snell \eqn{R^2}
#' @param shrinkage A pre-defined shrinkage factor. Default to 0.9
#' @param r2_nagelkerke An adjustment for Cox-Snell \eqn{R^2}. Default to 0.15.
#' @param auc A numeric vector. Pairwise C statistics
#' @param prev A numeric vector. Outcome proportions in categories relative to reference
#' @param sigma A double. Tolerate rate between adjusted \eqn{R^2} and apparent \eqn{R^2}. Default to 0.05
#' @param margin_error A double. A divergence for margin error specified by
#' Chi-square distribution's confidence level. Default to 0.05
#'
#' @return A \code{pmsample} class object
#' * \strong{final}: A single row data frame representing the selected criterion for
#' the minimum sample size for developing a new model
#' * \strong{criteria}: A data frame containing all the sample sizes
#' based on three criteria. Divergences in sample sizes due to number rounding
#' * sigma: A single numeric representing the tolerate rate between
#' adjusted \eqn{R^2} and apparent \eqn{R^2}
#' * margin_error: A divergence for margin error specified by
#' Chi-square distribution's confidence level
#' @export
#'
#' @examples
#' # when only C statistics is reported in a binary model
#' set.seed(1234)
#' pmsamplesize(Q = 30, k = 2, r2_nagelkerke = 0.15, auc = 0.81, prev = 0.77)
pmsamplesize <- function(Q, k, p, adjust = FALSE,
                         r2_cs_app = NULL, r2_cs_adj = NULL,
                         shrinkage = 0.9, r2_nagelkerke = 0.15,
                         auc = NULL, prev = NULL, sigma = 0.05,
                         margin_error = 0.05) {
  # constraints on k
  # inadequate information error msg
  err_msg <- "Inadequate information to derive `r2_cs_app`, try to supply a `r2_cs_app` or `r2_cs_adj` value"
  k <- ifelse(missing(k), length(p), k)
  if(missing(p)) {
    if(!is.null(r2_cs_app) | !is.null(r2_cs_adj)) {
      available_crit2 <- FALSE
    } else if(!is.null(auc) & !is.null(prev)) {
      available_crit2 <- FALSE
      if(k == 2) {
        r2_cs_app <- round(approximate_R2(k = k, auc = auc, prev = prev)[["R2.coxsnell"]], 2)
      } else {
        r2_cs_app <- mapply(approximate_R2, k = k, auc = auc, prev = prev)
      }
    }
  }
  stopifnot("k is not at least equal or greater than 2" = k >= 2)

  # k = 2
  if(k == 2) {

    ## constraints on p
    if(!missing(p)) is_prop <- all(p > 0 & p < 1)

    ## deriving r2_cs_adj: E, n information requires to be available
    if(is.null(r2_cs_adj)) {
      if(is.null(r2_cs_app)) {
        if(exists("is_prop")) {
          if(!is_prop) {
            ## when p are events' count
            ## the second value in the vector is the count of outcome
            n <- sum(p); E <- p[-1]; LR_null <- E*log(E/n) + (n-E)*log(1-E/n)
            # r2_cs_app = r2_nagelkerke * max(r2_cs_app)
            r2_cs_app <- r2_nagelkerke * round((1 - exp(2*LR_null/n)), 2)
            # if r2_cs_app is available prior
            # r2_cs_adj = s * r2_cs_app
            r2_cs_adj <- shrinkage * r2_cs_app
          } else stop(err_msg, call. = FALSE)
        }
      } else {
        r2_cs_adj <- shrinkage * r2_cs_app
      }
    } else if(!adjust) {
      r2_cs_app <- r2_cs_adj
    }

    ## whether r2_cs_adj = r2_cs_app?
    r2_cs_adj <- ifelse(adjust, r2_cs_adj, r2_cs_app)

    ### criterion 1: baseline on shrinkage = 0.9
    crit1 <- ceiling(Q / ((shrinkage - 1) * log(1 - r2_cs_adj/shrinkage)))
    ### criterion 2: baseline on adjusted shrinkage
    ### by absolute difference in proportions of variance explained
    if(!exists("available_crit2")) {
      shrinkage_crit2 <- round(r2_cs_adj / (r2_cs_adj + sigma*round(1 - exp(2*LR_null/n), 2)), 3)
      crit2 <- ceiling(Q / ((shrinkage_crit2 - 1) * log(1 - r2_cs_adj/shrinkage_crit2)))
    } else {
      warning("Inadequate information to calculate sample size for criterion 2")
      crit2 <- NA; shrinkage_crit2 <- NA
    }
    ### criterion 3: baseline on overall risk restrainted by margin error
    if(!is.null(prev)) {
      crit3 <- ceiling((1.96/margin_error)^2*prev*(1-prev))
    } else if(!missing(p)) {
      crit3 <- ceiling((1.96/margin_error)^2*(E/n)*(1-E/n))
    } else {
      warning("Inadequate information to calculate sample size for criterion 3")
      crit3 <- NA
    }

    res <- data.frame(
      sample_size = c(crit1, crit2, crit3),
      shrinkage = c(shrinkage, shrinkage_crit2, shrinkage),
      parameter = Q,
      CS_Rsq = r2_cs_adj,
      Max_Rsq = ifelse(missing(p), NA, round(1 - exp(2*LR_null/n), 2)),
      Nag_Rsq = r2_nagelkerke
    )

    if(missing(p)) {
      if(!is.null(prev)) {
        res$EPP <- round(res$sample_size * prev/Q, 2)
      } else {
        res$EPP <- NA
      }
    } else {
      res$EPP <- round((res$sample_size*E/n)/Q, 2)
    }

    res[4, ] <- res[which.max(res$sample_size), ]
    rownames(res) <- c(paste0("Criteria ", 1:3), "Final")
    if(!missing(p)) prevalance <- E/n
    if(!is.null(prev)) prevalance <- prev
    if(missing(p) & is.null(prev)) prevalance <- NA

    x <- list(final = res[4, ],
              criteria = res[complete.cases(res$sample_size), ],
              sigma = sigma,
              margin_error = margin_error,
              prev = prevalance)
    class(x) <- c("pmsample")
    return(x)
  }

  # k > 2
  if(k != 2) {
    prop <- all(p > 0 & p < 1)
    if(prop) {
      pk <- p; pkr <- margin.table(combn(p, 2), 2)
    } else {
      pk <- p / sum(p); pkr <- margin.table(combn(p, 2), 2) / sum(p)
    }

    if(!is.null(r2_cs_adj)) {
      if(length(r2_cs_adj) != ncol(combn(pk, 2))) {
        stop(sprintf('"r2_cs_adj" expects length of %g, but got %g',
                     ncol(combn(p, 2)), length(r2_cs_adj)))
      }
      r2_cs_adj_kr <- r2_cs_adj
    } else if(!is.null(r2_cs_app)) {
      if(length(r2_cs_app) != ncol(combn(pk, 2))) {
        stop(sprintf('"r2_cs_app" expects length of %g, but got %g',
                     ncol(combn(p, 2)), length(r2_cs_app)))
      }
      if(!adjust) {
        r2_cs_adj_kr <- r2_cs_app
      } else {
        r2_cs_adj_kr <- r2_nagelkerke * r2_cs_app
      }
    } else if(!is.null(auc)) {
      # throws an error if auc or prev is not the length of ncol(combn(p, 2))
      if(length(auc) != ncol(combn(pk, 2))) {
        stop(sprintf('"auc" expects length of %g, but got %g',
                     ncol(combn(p, 2)), length(auc)))
      }
      if(is.null(prev)) {
        prev <- pk[Reduce(c, sapply(seq_along(pk), function(i) seq_along(pk)[-seq_len(i)]))] / pkr
      } else if(length(prev) != ncol(combn(pk, 2))) {
        stop(sprintf('"prev" expects length of %g, but got %g',
                     ncol(combn(p, 2)), length(prev)))
      }
      r2_cs_adj_kr <- mapply(approximate_R2, k = k, auc = auc, prev = prev)
    } else {
      if(is.null(prev)) {
        prev <- pk[Reduce(c, sapply(seq_along(pk), function(i) seq_along(pk)[-seq_len(i)]))] / pkr
      } else if(length(prev) != ncol(combn(pk, 2))) {
        stop(sprintf('"prev" expects length of %g, but got %g',
                     ncol(combn(p, 2)), length(prev)))
      }
      # if neither pseudo-R2 or (pairwise) C-statistics are available a priori
      # estimates are typically larger than that from simulation
      r2_cs_app <- 1-((prev^prev)*((1-prev)^(1-prev)))^2
      r2_cs_adj_kr <- r2_nagelkerke * r2_cs_app
    }

    crit1 <- ceiling(max((Q / ((shrinkage - 1)*log(1-r2_cs_adj_kr/shrinkage))) / pkr))

    r2_cs_app <- 1 - prod(pk^pk)^2; r2_cs_adj <- r2_nagelkerke * r2_cs_app
    shrinkage_crit2 <- r2_cs_adj / (r2_cs_adj + sigma * r2_cs_app)
    crit2 <- ceiling(((k-1)*Q) / (shrinkage_crit2-1) / log(1-r2_cs_adj-sigma*r2_cs_app))

    crit3 <- ceiling(max(qchisq(0.05/k, 1, lower.tail = FALSE)*pk*(1-pk) / sigma^2))

    res <- data.frame(
      sample_size = c(crit1, crit2, crit3),
      shrinkage = c(shrinkage, shrinkage_crit2, shrinkage),
      parameter = Q,
      CS_Rsq = r2_cs_adj,
      Max_Rsq = r2_cs_app,
      Nag_Rsq = r2_nagelkerke
    )

    res[4, ] <- res[which.max(res$sample_size), ]
    rownames(res) <- c(paste0("Criteria ", 1:3), "Final")

    x <- list(final = res[4, ],
              criteria = res[complete.cases(res$sample_size), ],
              sigma = sigma,
              margin_error = margin_error)
    class(x) <- c("pmsample_multi")
    return(x)
  }
}
