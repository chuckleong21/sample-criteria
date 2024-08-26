#' Sample Criteria Extension for Multinomial Regression
#'
#' @param k An integer. The number of levels in response variable
#' @param Q An integer. The number of candidate predictors. See Details
#' @param props A numeric vector. The proportions of each events
#' @param n_events A integer vector. The number of each events
#' @param shrinkage Default is 0.9
#' @param r2_cs A vector
#' @param r2_nagelkerke
#' @param simulation
#' @param sim_params A list of pairwise C statistics and prevalences
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
pmsamplesize <- function(k,
                         Q = NULL,
                         props,
                         n_events,
                         shrinkage = 0.9,
                         r2_cs = NULL,
                         r2_nagelkerke = 0.15,
                         simulation = c("glm", "lrm"),
                         sim_params = NULL,
                         sigma = 0.05
) {


  # Calculating Q -----------------------------------------------------------
  # Since Q is available through fractional polynomial method and this
  # is modeled using the {mrp} package which requires empirical data
  # that is not currently available
  # For this work Q=17 will be used though personally I don't that's true
  # for different k
  Q <- Q %||% 17


  # Calculating p[k] & p[k,r] -----------------------------------------------

  if(missing(props) & missing(n_events)) {
    stop("One of `props` and `n_events` must be defined")
  } else if(!missing(props) & !missing(n_events)) {
    stop("Only one of `props` or `n_events` is used, not both")
  } else if(!missing(props) & missing(n_events)) {
    if(length(props) != k + ncol(combn(k, 2))) {
      stop(sprintf("Expect %d `props` values, but got %d", k + ncol(combn(k, 2)), length(props)))
    }
    n_events <- NULL
  } else {
    stopifnot(length(n_events) == k)
    props <- props_by_events(n_events)
  }


  # Calculating R-square[CS_app] -------------------------------------------
  if(is.null(names(props))) {
    kk <- Reduce(c, sapply(seq_len(k-1), function(i) seq(k)[-seq_len(i)]))
    rr <- rep(seq_len(k)[-k], rev(seq_len(k-1)))
    names(props) <- c(seq_len(k),
                  sprintf("%d_%d", kk, rr))
  }
  pk <- props[which(!grepl("_", names(props)))]
  r2_csapp <- 1 - (prod(pk^pk)^2)

  # Calculating R-square[CS_adj] --------------------------------------------
  r2_csadj <- r2_nagelkerke * r2_csapp

  # Calculating R-square[CS_adj,k,r] ----------------------------------------

  # if pairwise C-statistics are reported
  if(simulation %in% c("glm", "lrm")) {
    check_kwargs(sim_params)
    res <- mapply(approximate_R2,
                  method = simulation,
                  auc = sim_params$c_stats,
                  prev = sim_params$prev)

    if(identical("glm", simulation)) {
      names(res) <- NULL
    } else {
      res <- Reduce(c, res["R2.coxsnell", ])
    }
    r2_csadjkr <- res
  } else {
    message("`simulation` is either 'glm' or 'lrm'.")
    warning("argument `sim_args` is ignored for `simulation` argument is invalid")
    sim_args <- NULL
  }


  # criterion 1 -------------------------------------------------------------

  pkr <- props[which(grepl("_", names(props)))]

  if(!is.null(n_events)) {
    crit1 <- ceiling(max(Q / ((shrinkage-1)*log(1-r2_csadjkr/shrinkage)) * 1/pkr))
  } else {
    crit1 <- ceiling(max(Q / ((shrinkage-1)*log(1-r2_csadjkr/shrinkage)) * 1/pkr))
  }

# criterion 2 -------------------------------------------------------------

  crit2 <- ceiling((k-1) * Q / ((r2_csadj/(r2_csadj + sigma*r2_csapp)-1)*log(1-r2_csadj-sigma*r2_csapp)))

# criterion 3 -------------------------------------------------------------

  crit3 <- ceiling(max(qchisq(0.05/k, 1, lower.tail = FALSE)*pk*(1-pk) / sigma^2))

  max(c(crit1, crit2, crit3))
}
