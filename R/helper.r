phi_pairs <- function(n_events) {
  K <- length(n_events)
  x <- Reduce(c, sapply(seq_len(K-1), function(i) n_events[-seq_len(i)]))
  y <- rep(n_events[-K], c(rev(seq_len(K-1))))
  x / (x + y)
}

props_by_events <- function(x) {
  p <- c(x / sum(x),
         colSums(combn(x, 2)) / sum(x))
  k <- Reduce(c, sapply(seq_len(length(x)-1), function(i) seq(length(x))[-seq_len(i)]))
  r <- rep(seq_along(x)[-length(x)], rev(seq_len(length(x)-1)))
  names(p) <- c(seq_along(x),
                sprintf("%d_%d", k, r))
  list(events = x, p = p)
}

crit <- function(Q, shrinkage = 0.9,
                 r2_cs = NULL,
                 max_r2_csapp = NULL,
                 r2_csadj = NULL,
                 prev,
                 sigma = 0.05,
                 margin_error = 0.05,
                 ...) {
  x <- list(...)
  # criterion 1 -------------------------------------------------------------

  if("r2_csadjkr" %in% names(x)) {
    crit1 <- ceiling(max(Q / ((shrinkage-1)*log(1-x$r2_csadjkr/shrinkage)) * 1/x$pkr))
  } else {
    crit1 <- ceiling(Q/((shrinkage - 1)*log(1-r2_csadj/shrinkage)))
  }

  # criterion 2 -------------------------------------------------------------

  if("r2_csadjkr" %in% names(x)) {
    crit2 <- ceiling((x$k-1) * Q / ((r2_csadj/(r2_csadj + sigma*max_r2_csapp)-1)*log(1-r2_csadj-sigma*max_r2_csapp)))
  } else {
    crit2 <- ceiling(Q/((round(r2_csadj / (r2_csadj + sigma * max_r2_csapp), 3) - 1)*log(1-r2_csadj/(round(r2_csadj / (r2_csadj + sigma * max_r2_csapp), 3)))))
  }
  # criterion 3 -------------------------------------------------------------

  if("pkr" %in% names(x)) {
    crit3 <- ceiling(max(qchisq(0.05/x$k, 1, lower.tail = FALSE)*x$pk*(1-x$pk) / sigma^2))
  } else {
    crit3 <- ceiling((1.96/margin_error)^2*prev*(1-prev))
  }

  c(crit1, crit2, crit3)
}
