test_that("Only C statistics & outcome proportion are reported", {
  set.seed(1234)
  s <- suppressWarnings(pmsamplesize(Q = 30, k = 2, r2_nagelkerke = 0.15,
                    auc = 0.81, prev = 0.77))
  expect_equal(s$final$sample_size, 1130)
})

test_that("When R2_CS is not reported", {
  s <- pmsamplesize(Q = 24, k = 2, p =c(138-24, 24), r2_nagelkerke = 0.48, auc = 0.91)
  expect_equal(s$final$sample_size, 668)
  expect_no_error(pmsamplesize(Q = 24, k = 2, p =c(138-24, 24), r2_nagelkerke = 0.48, auc = 0.91))
})

test_that("Inadequate information error", {
  # sample size can't be calculated since neither E, n are known
  expect_error(pmsamplesize(Q = 24, p = c(0.826, 0.174), r2_nagelkerke = 0.48))
  v <- suppressWarnings(
    pmsamplesize(Q = 24, k = 2, r2_cs_app = 0.288, prev = 0.174, r2_nagelkerke = 0.48)
  )
  vv <- suppressWarnings(
    pmsamplesize(Q = 24, k = 2, r2_cs_app = 0.288, r2_nagelkerke = 0.48)
  )
  # Only criterion 1 is calculated
  expect_equal(rownames(vv$criteria)[-nrow(vv$criteria)], "Criteria 1")
  # Only criterion 1 & 3 are calculated
  expect_equal(rownames(v$criteria)[-nrow(v$criteria)], c("Criteria 1", "Criteria 3"))
})

test_that("Only a certain criteria is calculated", {
  s <- suppressWarnings(pmsamplesize(Q = 24, k = 2, r2_cs_app = 0.288, r2_nagelkerke = 0.48))
  expect_equal(s$final$sample_size, 623)
})

test_that("Results are different based on the information type of p", {
  set.seed(101)
  a <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(0.729, 0.053, 0.05, 0.133, 0.034),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9,
                    auc = c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82))
  set.seed(101)
  b <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(2557, 186, 176, 467, 120),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9,
                    auc = c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82))
  set.seed(101)
  x <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(2557, 186, 176, 467, 120),
                    r2_cs_adj = c(0.116, 0.179, 0.497, 0.170, 0.185, 0.499, 0.374, 0.328, 0.129, 0.210),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9)
  set.seed(101)
  y <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(2557, 186, 176, 467, 120),
                    r2_cs_app = c(0.391, 0.38, 0.577, 0.306, 0.75, 0.697, 0.738, 0.691, 0.741, 0.637),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9)
  expect_equal(a$final$sample_size, 13135)
  expect_equal(b$final$sample_size, 13063)
  expect_equal(x$final$sample_size, 13016)
  expect_equal(y$final$sample_size, 15276)
})


