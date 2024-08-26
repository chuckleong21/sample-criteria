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
  p
}

check_kwargs <- function(kwargs) {
  arg_match <- lapply(names(kwargs), function(x) try(match.arg(x, c("c_stats", "prev")), silent = TRUE))
  stopifnot("One or more arguments are not `c_stats`, `prev` when simulation is FALSE" = !all(sapply(arg_match, function(x) inherits(x, "try-error"))))
  invisible()
}
