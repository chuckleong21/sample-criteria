set.seed(101)
c_stats <- c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82)
events <- c(2557, 186, 176, 467, 120)
prev <- phi_pairs(events)

test_that("A legal k should be scalar and an integer that is at least 2", {
  set.seed(1234)
  expect_error(approximate_R2(k = 1, auc = c_stats, prev = prev))
  expect_warning(approximate_R2(k = 2.5, auc = 0.81, prev = 0.77))
  expect_error(approximate_R2(k = 1, auc= 0.81, prev = 0.77))
})

test_that("Approximate R2 based on binary outcome", {
  set.seed(1234)
  r2 <- lapply(approximate_R2(k = 2, auc = 0.81, prev = 0.77), function(x) round(x, 3))
  expect_equal(r2$R2.nagelkerke, round(0.3183689, 3))
  expect_equal(r2$R2.coxsnell, round(0.2100957, 3))
})

test_that("Approximate R2 based on multinomial outcome", {
  set.seed(101)
  r2 <- round(mapply(approximate_R2, k = 5, auc = c_stats, prev = prev, USE.NAMES = FALSE), 3)
  expected <- round(c(0.1161076, 0.1793267, 0.4967071, 0.1695518, 0.1851988, 0.4993293, 0.3738221, 0.3280985, 0.1285719, 0.2104207), 3)
  expect_equal(r2, expected)
})
