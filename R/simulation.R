#' R-sqaure approximation
#'
#' @description
#' \eqn{R^2} approximation based on Riley's approach when \eqn{R^2_{CS}} is not
#' reported.
#'
#'
#' @importFrom rms lrm
#' @param k An scalar integer. The number of outcomes, \eqn{k\geq{2}}
#' @param auc A double. Reported pairwise C-statistics
#' @param prev A double. The outcome  reporting the existing prediction model
#' @param n An integer. Simulation sample size. Default to \emph{1e6}
#'
#' @return
#' * For binary outcome, approximate_R2() bases on \insertCite{Riley2021-qe}{sample.criteria}
#' using \code{rms} package. It returns a list containing the \eqn{R^2_{CS}} and \eqn{R^2_{Nagelkerke}}
#'
#' * While for multinomial models, A Bernoulli simulation is used, based on
#' \insertCite{Pate2023-yh}{sample.criteria}
#' @export
#'
#' @examples
#' sim <- approximate_R2(k = 2, auc = 0.81, prev = 0.77)
#' sim[["R2.coxsnell"]] # coxsnell R2
#' sim[["R2.nagelkerke"]] # nagerlkerke R2
approximate_R2 <- function(k, auc, prev, n = 1e6) {
  # Make sure k is an integer
  stopifnot("k must be scalar" = length(k) == 1)
  if(k %% 1 != 0) warning(paste0("k = ", deparse(substitute(k)), ", Coerce k = ", as.integer(k)))
  k <- as.integer(k)

  # binary outcome model
  if(k == 2) {
    # define mu as a function of the C statistic
    mu <- sqrt(2) * qnorm(auc)

    # simulate large sample linear prediction based on two normals
    # for non-eventsN(0, 1), events and N(mu, 1)

    LP <- c(rnorm(prev * n, mean = 0, sd = 1), rnorm((1 - prev) * n, mean = mu, sd = 1))
    y <- c(rep(0, prev * n), rep(1, (1 - prev) * n))

    # Fit a logistic regression with LP as covariate;
    # this is essentially a calibration model, and the intercept and
    # slope estimate will ensure the outcome proportion is accounted
    # for, without changing C statistic

    fit <- lrm(y ~ LP)

    max_R2 <- function(prev) {
      1 - (prev^prev * (1 - prev)^(1 - prev))^2
    }

    return(list(
      R2.nagelkerke = as.numeric(fit$stats["R2"]),
      R2.coxsnell = as.numeric(fit$stats["R2"]) * max_R2(prev)
    ))
  } else if(k > 2) {
    # Create an empty dataset
    out <- data.frame(matrix(ncol = 2, nrow = n))
    colnames(out) <- c("Y", "LP")

    # Create the outcome variable
    Y_vec <- rbinom(n, 1, prev)

    # Create the vector of mean values for the linear predictor data generation
    Y_vec_mu <- Y_vec*sqrt(2)*qnorm(auc, 0, 1)

    # Generate the linear predictor
    LP_vec <- rnorm(n, Y_vec_mu, 1)

    # Assign these into an output dataset
    out$Y <- as.integer(Y_vec)
    out$LP <- LP_vec

    # general model
    model <- glm(Y ~ LP_vec, family = binomial(link = "logit"), data = out)
    # null model
    model_null <- glm(Y ~ 1, family = binomial(link = "logit"), data = out)

    # Calculate likelihood ratio statistics
    LR <- as.numeric(-2*(logLik(model_null) - logLik(model)))

    # Calculate r2cs_app
    r2cs_app <- 1 - exp(-LR/n)
    r2cs_app
  } else {
    stop("k must be at least 2")
  }
}

