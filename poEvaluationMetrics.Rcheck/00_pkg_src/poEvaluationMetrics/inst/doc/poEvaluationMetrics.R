## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----helpers------------------------------------------------------------------
library(poEvaluationMetrics)

actual <- c(1, 1, 0, 0)
pred_bin <- c(1, 0, 0, 1)

Fbp(actual, pred_bin)
omission(actual, pred_bin)
orss(actual, pred_bin)
brier_score(c(0, 1, 1, 0), c(0.2, 0.7, 0.9, 0.1))

## ----indexCalculation, message = FALSE----------------------------------------
library(terra)
library(sf)

set.seed(1)
r <- rast(
  nrows = 80,
  ncols = 80,
  xmin = 0,
  xmax = 1,
  ymin = 0,
  ymax = 1,
  crs = "local"
)
values(r) <- runif(ncell(r))

inputDF <- data.frame(
  predicted = c(0.85, 0.72, 0.25, 0.41, 0.05, 0.91),
  observed = c(1L, 1L, 0L, 0L, 0L, 1L)
)

round(indexCalculation(inputDF, prediction = r), 4)

## ----calculateMetrics---------------------------------------------------------
pres_coords <- data.frame(x = c(0.2, 0.35, 0.6), y = c(0.25, 0.7, 0.45))
bg_coords   <- data.frame(x = c(0.1, 0.5, 0.85), y = c(0.55, 0.15, 0.75))

presence <- st_as_sf(pres_coords, coords = c("x", "y"), crs = crs(r))
absence  <- st_as_sf(bg_coords, coords = c("x", "y"), crs = crs(r))

round(calculateMetrics(r, presence, absence), 4)

## ----performanceEstimation, message = FALSE-----------------------------------
env <- c(r, r * 0.8 + 0.1)
names(env) <- c("v1", "v2")

pe_result <- performanceEstimation(
  prediction = r,
  presence = presence,
  absence = absence,
  background = TRUE,
  aa = FALSE,
  environmentalVariables = env,
  noPointsTesting = 20,
  replicates = 2
)

pe_result[, 1:8]

