# source("R/helper.r")
# source("R/simulation.R")

K <- 5
events <- c(2557, 186, 176, 467, 120)
c_stats <- c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82)
prev <- phi_pairs(events)

test_that("Test Arguments: Only one of props or events should be supplied.", {
  p <- events / sum(events)
  expect_error(pmsamplesize(k = K), "One of `props` and `n_events` must be defined")
  expect_error(
    pmsamplesize(k = K,
                 props = props_by_events(events),
                 n_events = events,
                 simulation = "glm",
                 sim_params = list(c_stats = c_stats, prev = prev)),
    "Only one of `props` or `n_events` is used, not both"
  )
  expect_error(pmsamplesize(k = K, props = p),
               "Expect 15 `props` values, but got 5")
})

test_that("Either specification should return the same result", {
  set.seed(101)
  p <- props_by_events(events)
  sample_by_props <- pmsamplesize(k = K, props = p,
                                  simulation = "glm",
                                  sim_params = list(c_stats = c_stats, prev = prev))
  set.seed(101)
  sample_by_events <- pmsamplesize(k = K, n_events = events,
                                   simulation = "glm",
                                   sim_params = list(c_stats = c_stats, prev = prev))
  expect_identical(sample_by_props,sample_by_events)
})


