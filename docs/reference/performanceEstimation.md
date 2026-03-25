# SDM performance evaluation

Evaluates the performance of species distribution models (SDMs). The
evaluation supports three data types: presence-absence (PA),
presence-artificial-absence (PAA), and presence-background (PBG).
Metrics include predictive accuracy (AUC, COR), predictive error
measures (MAE, BIAS), across cross-validation replicates.

## Usage

``` r
performanceEstimation(
  prediction,
  presence,
  absence = FALSE,
  background = TRUE,
  aa = TRUE,
  environmentalVariables = NA,
  noPointsTesting = NA,
  replicates = 100
)
```

## Arguments

- prediction:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the prediction map.

- presence:

  **Required.** An `sf` object of known presence locations.

- absence:

  Optional. An `sf` object of known absence locations, or `FALSE` to
  skip PA metrics.

- background:

  Optional. An `sf` object of background points, or logical `TRUE` to
  auto-generate, or `FALSE` to skip PBG.

- aa:

  Optional. An `sf` object of artificial absence points, or logical
  `TRUE` to derive via AOA, or `FALSE` to skip PAA.

- environmentalVariables:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of covariates when `background` or `aa` is `TRUE`, or `sf` points are
  not supplied for those; may be `NA` only if both are `FALSE` or only
  fixed `sf` points are used.

- noPointsTesting:

  Integer. Number of background or artificial absence points to generate
  when sampling.

- replicates:

  Integer. Number of replicate draws for PBG/PAA; column-wise means
  across replicates are returned for those scenarios.

## Value

A `data.frame` with one row per computed scenario. The first column is
`scenario` (`"PA"`, `"PAA"`, or `"PBG"`). Remaining columns match
[`calculateMetrics()`](https://envima.github.io/poEvaluationMetrics/reference/calculateMetrics.md)
/
[`indexCalculation()`](https://envima.github.io/poEvaluationMetrics/reference/indexCalculation.md).

## Details

Metrics combine discrimination (e.g. AUC, COR), probability calibration
errors (MAE, BIAS), and threshold-based skill from
[`mecofun::evalSDM()`](https://rdrr.io/pkg/mecofun/man/evalSDM.html) and
helper indices
([`Fbp()`](https://envima.github.io/poEvaluationMetrics/reference/Fbp.md),
[`sedi()`](https://envima.github.io/poEvaluationMetrics/reference/sedi.md),
[`orss()`](https://envima.github.io/poEvaluationMetrics/reference/orss.md),
[`omission()`](https://envima.github.io/poEvaluationMetrics/reference/omission.md)).
PAA uses CAST `aoa()` to mask environmentally dissimilar areas before
sampling artificial absences.

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- performanceEstimation(
    prediction = prediction_raster,
    presence = presence_points,
    background = TRUE,
    environmentalVariables = env_rasters,
    replicates = 50
  )
  result[result$scenario == "PBG", ]
} # }
```
