skip_integration <- function() {
  pkgs <- c(
    "terra", "sf", "mecofun", "PresenceAbsence", "prg", "Metrics",
    "CAST", "predicts", "mgcv", "dplyr"
  )
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      testthat::skip(sprintf("package '%s' not installed", p))
    }
  }
}

test_that("indexCalculation runs on small synthetic data", {
  skip_integration()

  r <- terra::rast(
    nrows = 60, ncols = 60,
    xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    crs = "local"
  )
  terra::values(r) <- stats::runif(terra::ncell(r))

  inputDF <- data.frame(
    predicted = c(0.85, 0.72, 0.25, 0.41, 0.05, 0.91),
    observed = c(1L, 1L, 0L, 0L, 0L, 1L)
  )

  out <- indexCalculation(inputDF, prediction = r)
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true("AUC" %in% names(out))
  testthat::expect_true(nrow(out) == 1L)
})

test_that("calculateMetrics runs with sf presence and absence", {
  skip_integration()

  r <- terra::rast(
    nrows = 60, ncols = 60,
    xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    crs = "local"
  )
  terra::values(r) <- stats::runif(terra::ncell(r))

  pres <- sf::st_as_sf(
    data.frame(x = c(0.2, 0.35, 0.6), y = c(0.25, 0.7, 0.45)),
    coords = c("x", "y"),
    crs = terra::crs(r)
  )
  absn <- sf::st_as_sf(
    data.frame(x = c(0.1, 0.5, 0.85), y = c(0.55, 0.15, 0.75)),
    coords = c("x", "y"),
    crs = terra::crs(r)
  )

  out <- calculateMetrics(r, pres, absn)
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true(nrow(out) == 1L)
})

test_that("performanceEstimation returns data.frame with scenario column", {
  skip_integration()

  r <- terra::rast(
    nrows = 50, ncols = 50,
    xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    crs = "local"
  )
  terra::values(r) <- stats::runif(terra::ncell(r))
  env <- c(r, r * 0.8 + 0.1)
  names(env) <- c("v1", "v2")

  pres <- sf::st_as_sf(
    data.frame(x = c(0.3, 0.55), y = c(0.4, 0.6)),
    coords = c("x", "y"),
    crs = terra::crs(r)
  )
  absn <- sf::st_as_sf(
    data.frame(x = c(0.15, 0.8), y = c(0.7, 0.25)),
    coords = c("x", "y"),
    crs = terra::crs(r)
  )

  out <- performanceEstimation(
    prediction = r,
    presence = pres,
    absence = absn,
    background = TRUE,
    aa = FALSE,
    environmentalVariables = env,
    noPointsTesting = 15,
    replicates = 2
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true("scenario" %in% names(out))
  testthat::expect_true(all(out$scenario %in% c("PA", "PBG", "PAA")))
  testthat::expect_true(nrow(out) >= 1L)
})
