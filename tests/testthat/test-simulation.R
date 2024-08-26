# source("R/helper.r")

set.seed(101)
c_stats <- c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82)
events <- c(2557, 186, 176, 467, 120)
prev <- phi_pairs(events)
r2cs_adjkr_glm <- mapply(approximate_R2, method = "glm", auc = c_stats, prev = prev)
r2cs_adjkr_lrm <- mapply(approximate_R2, method = "lrm", auc = c_stats, prev = prev)
names(r2cs_adjkr_glm) <- NULL
r2cs_adjkr_lrm <- Reduce(c, r2cs_adjkr_lrm["R2.coxsnell", ])


test_that("Method test: glm", {
  r2cs_adjkr_expected <- c(0.1161076, 0.1793267, 0.4967071, 0.1695518, 0.1851988, 0.4993293, 0.3738221, 0.3280985, 0.1285719, 0.2104207)
  expect_equal(round(r2cs_adjkr_glm, 3), round(r2cs_adjkr_expected, 3))
})

test_that("Method test: divergence between two methods all else equal", {
  expect_failure(expect_equal(r2cs_adjkr_glm, r2cs_adjkr_lrm))
})

