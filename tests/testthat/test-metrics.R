test_that("Fbp matches contingency definition", {
  actual <- c(1, 1, 0, 0)
  pred <- c(1, 0, 1, 0)
  tp <- 1L
  fp <- 1L
  fn <- 1L
  expect_equal(Fbp(actual, pred), (2 * tp) / (tp + fn + fp))
})

test_that("omission is FN / (TP + FN)", {
  actual <- c(1, 1, 0)
  pred <- c(0, 1, 0)
  expect_equal(omission(actual, pred), 1 / 2)
})

test_that("brier_score validates length", {
  expect_error(brier_score(c(0, 1), c(0.1)), "Vectors must be same length")
  expect_equal(brier_score(c(0, 1), c(0.2, 0.8)), mean(c(0.2^2, 0.2^2)))
})

test_that("orss is odds ratio skill score", {
  actual <- c(1, 1, 0, 0)
  pred <- c(1, 0, 0, 1)
  tp <- 1L
  tn <- 1L
  fp <- 1L
  fn <- 1L
  expect_equal(orss(actual, pred), (tp * tn - fp * fn) / (tp * tn + fp * fn))
})
