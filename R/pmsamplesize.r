#' @exportS3Method base::print
print.pmsample <- function(x, ...) {
  message(paste("NB: Assuming", x$sigma, "acceptable difference in apparent & adjusted R-squared"))
  message(paste("NB: Assuming", x$margin_error, "margin of error in estimation of intercept"))
  message(paste("NB: Events per Predictor Parameter assumes prevalance", round(x$prev, 3), "\n"))
  print(x$criteria)
  cat(sprintf("Minimum sample size required for new model development based on user inputs = %.0f,
 with %.0f events (assuming an outcome prevalence = %.3g) and an EPP = %.3g",
              x$final$sample_size,
              x$final$sample_size * x$prev,
              x$prev, x$final$EPP))
}

#' @exportS3Method base::print
print.pmsample_multi <- function(x, ...) {
  message(paste("NB: Assuming", x$sigma, "acceptable difference in apparent & adjusted R-squared"))
  message(paste("NB: Assuming", x$margin_error, "margin of error in estimation of intercept"))
  print(x$criteria)
}

pmsamplesize_binary <- function(Q,
                                p,
                                r2_cs,
                                shrinkage = 0.9,
                                r2_nagelkerke = 0.15,
                                auc,
                                prev,
                                sigma = 0.05,
                                margin_error = 0.05) {
  if(is.null(r2_cs) & missing(p)) {
    r2 <- approximate_R2(k = 2, auc = auc, prev = prev)
    max_r2_csapp <- r2[["R2.nagelkerke"]]
    r2_csadj <- r2[["R2.coxsnell"]]

    crits <- crit(Q = Q,
                  max_r2_csapp = max_r2_csapp,
                  r2_csadj = r2_csadj,
                  prev = prev)
  }

  if(is.null(r2_cs) & !missing(p)) {
    k <- 2; events <- p$events; pr <- p$p[-which(p$p == 1)]
    names(pr) <- NULL; prev <- ifelse(missing(prev), pr[k], prev)
    LR <- round(events[k]*log(pr[k])+(sum(events)-events[k])*log(1-pr[k]), 2)
    max_r2_csapp <- round(1-exp(2*LR/sum(events)), 2)
    r2_csapp <- r2_csadj <- r2_nagelkerke * max_r2_csapp

    crits <- crit(Q = Q,
                  max_r2_csapp = max_r2_csapp,
                  r2_csadj = r2_csadj,
                  prev = prev)
  }

  if(!is.null(r2_cs)) {
    max_r2_csapp <- r2_cs; r2_csadj <- r2_nagelkerke * max_r2_csapp
    crits <- crit(Q = Q,
                  max_r2_csapp = max_r2_csapp,
                  r2_csadj = r2_csadj,
                  prev = prev)
  }

  res <- data.frame(sample_size = crits,
                    shrinkage = c(shrinkage, rep(r2_csadj / (r2_csadj + sigma * max_r2_csapp), 2)),
                    parameter = Q,
                    CS_Rsq = r2_csadj,
                    Max_Rsq = max_r2_csapp,
                    Nag_Rsq = ifelse(missing(p), max_r2_csapp, r2_nagelkerke),
                    EPP = crits * prev/Q)
  res[4, ] <- res[which.max(res$sample_size), ]
  rownames(res) <- c(paste0("Criteria ", 1:3), "Final")

  x <- list(final = res[4, ],
            criteria = res,
            sigma = sigma,
            margin_error = margin_error,
            prevalance = prev)
  class(x) <- c("pmsample")
  x
}

pmsamplesize_multi <- function(Q,
                               k,
                               p,
                               r2_cs,
                               shrinkage = 0.9,
                               r2_nagelkerke = 0.15,
                               auc = NULL,
                               prev = NULL,
                               sigma = 0.05,
                               margin_error = 0.05) {
  events <- p$events; pr <- p$p
  if(is.null(names(pr))) {
    kk <- Reduce(c, sapply(seq_len(k-1), function(i) seq(k)[-seq_len(i)]))
    rr <- rep(seq_len(k)[-k], rev(seq_len(k-1)))
    names(pr) <- c(seq_len(k),
                   sprintf("%d_%d", kk, rr))
  }

  pk <- pr[which(!grepl("_", names(pr)))]
  pkr <- pr[-which(!grepl("_", names(pr)))]

  if(!is.null(prev)) {
    prev
  } else {
    prev <- phi_pairs(events)
  }
  max_r2_csapp <- 1 - (prod(pk^pk)^2)
  r2_csadj <- r2_nagelkerke * max_r2_csapp
  if(is.null(auc)) {
    phis <- phi_pairs(events)
    max_r2_csappkr <- 1 - (phis^phis*(1-phis)^(1-phis))^2
    r2_csadjkr <- r2_nagelkerke*max_r2csappkr
  }
  r2_csadjkr <- mapply(approximate_R2, k = k, auc = auc, prev = prev)
  crits <- crit(Q = Q,
                max_r2_csapp = max_r2_csapp,
                r2_csadj = r2_csadj,
                r2_csadjkr = r2_csadjkr,
                k = k,
                pk = pk,
                pkr = pkr)
  res <- data.frame(sample_size = crits,
                    shrinkage = c(shrinkage, rep(r2_csadj / (r2_csadj + sigma * max_r2_csapp), 2)),
                    parameter = Q,
                    CS_Rsq = r2_csadj,
                    Max_Rsq = max_r2_csapp,
                    Nag_Rsq = ifelse(missing(p), max_r2_csapp, r2_nagelkerke)
                    # EPP = crits * prev/Q
  )
  res[4, ] <- res[which.max(res$sample_size), ]
  rownames(res) <- c(paste0("Criteria ", 1:3), "Final")

  x <- list(final = res[4, ],
            criteria = res,
            sigma = sigma,
            margin_error = margin_error,
            prevalance = prev)
  class(x) <- "pmsample_multi"
  x
}

#' Sample Criteria Extension for Multinomial Regression
#'
#' @importFrom Rdpack reprompt
#'
#' @description
#' Calculate sample criteria extended for multinominal regression.
#'
#' @details
#' Sample Criteria for Clinical Prediction Model are based on the work of
#' \href{https://onlinelibrary.wiley.com/doi/full/10.1002/sim.7992}{Riley 2018} and
#' \href{https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.8806}{Calster 2020}
#'
#' @references
#' \insertCite{pmsamplesize}{sample.criteria}
#'
#' @param Q An integer. The number of candidate predictors. See Details
#' @param k An integer. The number of levels in response variable
#' @param p A numeric vector. Outcome proportion(s)
#' @param shrinkage A pre-defined shrinkage factor. Default to 0.9
#' @param r2_cs numeric. Cox-Snell \eqn{R^2}.
#' @param r2_nagelkerke An adjustment for Cox-Snell \eqn{R^2}. Default to 0.15.
#' @param auc numeric. Pairwise C statistics
#' @param prev numeric. Outcome proportions in categories relative to reference
#' @param sigma A double. tolerate rate between adjusted \eqn{R^2} and apparent \eqn{R^2}. Default to 0.05
#' @param margin_error A divergence for margin error specified by Chi-sqaure distribution's confidence level.
#' Deault to 0.05
#'
#' @return An integer
#' @export
#'
#' @examples
#' K <- 5
#' events <- c(2557, 186, 176, 467, 120)
#' c_stats <- c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82)
#' prev <- phi_pairs(events)
#' set.seed(1234)
#' pmsamplesize(Q = 30, k = 2, auc = 0.81, prev = 0.77)
pmsamplesize <- function(Q, k, p,
                         r2_cs = NULL,
                         auc = NULL,
                         prev,
                         shrinkage  = 0.9,
                         r2_nagelkerke = 0.15,
                         sigma = 0.05,
                         margin_error = 0.05) {
  if(k == 2) {
    pmsamplesize_binary(Q = Q, p = p, r2_cs = r2_cs,
                        shrinkage = shrinkage, r2_nagelkerke = r2_nagelkerke, auc = auc, prev = prev,
                        sigma = sigma, margin_error = margin_error)
  } else {
    pmsamplesize_multi(Q = Q, k = k, p = p, r2_cs = r2_cs,
                       shrinkage = shrinkage, r2_nagelkerke = r2_nagelkerke, auc = auc, prev = prev,
                       sigma = sigma, margin_error = margin_error)
  }
}
