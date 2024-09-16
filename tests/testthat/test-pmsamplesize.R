test_that("Valid k argument", {
  expect_error(pmsamplesize(Q = 30, k = 1, auc = 0.81, prev = 0.77))
})

test_that("Only C statistics & outcome proportion are reported", {
  set.seed(1234)
  s <- suppressWarnings(pmsamplesize(Q = 30, k = 2, r2_nagelkerke = 0.15,
                                     auc = 0.81, prev = 0.77))
  expect_equal(s$final$sample_size, 1130)
})

test_that("When R2_CS is not reported", {
  s <- pmsamplesize(Q = 24, k = 2, p = c(138-24, 24), r2_nagelkerke = 0.48, auc = 0.91)
  expect_equal(s$final$sample_size, 668)
  expect_no_error(pmsamplesize(Q = 24, k = 2, p = c(138-24, 24), r2_nagelkerke = 0.48, auc = 0.91))
})

test_that("Inadequate information error", {
  # sample size can't be calculated since neither E, n are known
  expect_error(pmsamplesize(Q = 24, p = c(0.826, 0.174), r2_nagelkerke = 0.48))

})

test_that("Only a certain criteria is calculated", {
  # EPP is unable to be calculated for no prevalance information
  s <- suppressWarnings(pmsamplesize(Q = 24, k = 2, r2_cs_app = 0.288, r2_nagelkerke = 0.48))
  # Only criterion 1 is calculated
  v <- suppressWarnings(
    pmsamplesize(Q = 24, k = 2, r2_cs_app = 0.288, prev = 0.174, r2_nagelkerke = 0.48)
  )
  vvv <- suppressWarnings(
    pmsamplesize(Q = 24, k = 2, r2_cs_adj = 0.288, r2_nagelkerke = 0.48)
  )
  # Only criterion 1 & 3 are calculated
  vv <- suppressWarnings(
    pmsamplesize(Q = 24, k = 2, r2_cs_app = 0.288, r2_nagelkerke = 0.48)
  )
  expect_equal(rownames(vvv$criteria)[-nrow(vvv$criteria)], "Criteria 1")
  expect_equal(rownames(vv$criteria)[-nrow(vv$criteria)], "Criteria 1")
  expect_equal(rownames(v$criteria)[-nrow(v$criteria)], c("Criteria 1", "Criteria 3"))
  expect_equal(s$final$sample_size, 623)
})

test_that("Valid amount of r2_cs, auc & prev", {
  expect_error(
    pmsamplesize(Q = 17,
                 k = 5,
                 p = c(2557, 186, 176, 467, 120),
                 r2_cs_app = c(0.391, 0.38, 0.306, 0.75, 0.697, 0.738, 0.691, 0.741, 0.637),
                 r2_nagelkerke = 0.15,
                 shrinkage = 0.9),
    "expects length"
  )
  expect_error(
    pmsamplesize(Q = 17,
                 k = 5,
                 p = c(2557, 186, 176, 467, 120),
                 r2_cs_adj = c(0.116, 0.179, 0.497, 0.170, 0.499, 0.374, 0.328, 0.129, 0.210),
                 r2_nagelkerke = 0.15,
                 shrinkage = 0.9),
    "expects length"
  )

  expect_error(
    pmsamplesize(Q = 17,
                 k = 5,
                 p = c(0.729, 0.053, 0.05, 0.133, 0.034),
                 r2_nagelkerke = 0.15,
                 shrinkage = 0.9,
                 auc = c(0.85, 0.92, 0.99, 0.95, 0.95, 0.87, 0.87, 0.71, 0.82)),
    "expects length"
  )
  expect_error(
    pmsamplesize(Q = 17,
                 k = 5,
                 p = c(0.729, 0.053, 0.05, 0.133, 0.034),
                 r2_nagelkerke = 0.15,
                 shrinkage = 0.9,
                 auc = c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82),
                 prev = c(0.064, 0.154, 0.045, 0.485, 0.715, 0.391, 0.727, 0.405, 0.204)),
    "expects length"
  )
})

test_that("Results are different based on the information type of p", {
  # when proportions & C-statistics are provided
  set.seed(101)
  a <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(0.729, 0.053, 0.05, 0.133, 0.034),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9,
                    auc = c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82))
  # when event outcomes & C-statistics are given
  set.seed(101)
  b <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(2557, 186, 176, 467, 120),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9,
                    auc = c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82))
  # when outcome events and adjusted cox-snell r2 are given
  set.seed(101)
  x <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(2557, 186, 176, 467, 120),
                    r2_cs_adj = c(0.116, 0.179, 0.497, 0.170, 0.185, 0.499, 0.374, 0.328, 0.129, 0.210),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9)
  # when outcome events and apparent cox-snell r2 are given
  y <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(2557, 186, 176, 467, 120),
                    adjust = TRUE, # use adjusted R2
                    r2_cs_app = c(0.391, 0.38, 0.577, 0.306, 0.75, 0.697, 0.738, 0.691, 0.741, 0.637),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9)
  # when only proportions are obtained
  z <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(0.729, 0.053, 0.05, 0.133, 0.034),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9)
  # when only outcome events are provided
  o <- pmsamplesize(Q = 17,
                    k = 5,
                    p = c(2557, 186, 176, 467, 120),
                    r2_nagelkerke = 0.15,
                    shrinkage = 0.9)
  expect_equal(a$final$sample_size, 13135)
  expect_equal(b$final$sample_size, 13063)
  expect_equal(x$final$sample_size, 13016)
  expect_equal(y$final$sample_size, 15276)
  expect_equal(z$final$sample_size, 15360)
  expect_equal(o$final$sample_size, 15280)
})


