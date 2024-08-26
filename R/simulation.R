#' Title
#'
#' @importFrom rms lrm
#' @param auc
#' @param prev
#' @param n
#'
#' @return
#' @export
#'
#' @examples
approximate_R2 <- function(method, auc, prev, n = 1e6) {
  switch (method,
    "glm" = {
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
    },
    "lrm" = {
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

      fit <- rms::lrm(y ~ LP)

      max_R2 <- function(prev) {
        1 - (prev^prev * (1 - prev)^(1 - prev))^2
      }

      return(list(
        R2.nagelkerke = as.numeric(fit$stats["R2"]),
        R2.coxsnell = as.numeric(fit$stats["R2"]) * max_R2(prev)
      ))
    }
  )
}

