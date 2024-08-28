events <- c(2557, 186, 176, 467, 120)
c_stats <- c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82)

test_that("Only C statistics & outcome proportion are reported", {
  set.seed(1234)
  expect_no_error(pmsamplesize(Q = 30, k = 2, auc = 0.81, prev = 0.77))
})

test_that("When R2_CS is not reported", {
  expect_no_error(pmsamplesize(Q = 24, k = 2, p = props_by_events(c(138-24, 24)), r2_nagelkerke = 0.48, auc = 0.91))
})

test_that("Simulation R2 when no r2_csadjk are reported", {
  set.seed(101)
  expect_no_error(pmsamplesize(Q = 17, k = 5, p = props_by_events(events), auc = c_stats, prev = NULL))
})


