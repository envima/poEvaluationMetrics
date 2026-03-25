# Getting started with poEvaluationMetrics

This package evaluates species distribution models (SDMs) with
discrimination metrics (AUC, correlation, PRG), calibration error (MAE,
bias), Boyce-type indices (`SBI_*`), and threshold-based skill scores
(TSS, Kappa, `Fbp`, SEDI, ORSS, omission, and related helpers).

## Standalone metric helpers

Many helpers take **binary** `actual` and `predicted` vectors (`0`/`1`):

``` r

library(poEvaluationMetrics)

actual <- c(1, 1, 0, 0)
pred_bin <- c(1, 0, 0, 1)

Fbp(actual, pred_bin)
#> [1] 0.6666667
omission(actual, pred_bin)
#> [1] 0.5
orss(actual, pred_bin)
#> [1] 0
brier_score(c(0, 1, 1, 0), c(0.2, 0.7, 0.9, 0.1))
#> [1] 0.0375
```

## From probabilities to a full metric table

[`indexCalculation()`](https://envima.github.io/poEvaluationMetrics/reference/indexCalculation.md)
expects a `data.frame` with columns `predicted` (numeric scores) and
`observed` (`0`/`1`), plus a `SpatRaster` of predictions (used for
Boyce/SBI sampling and related steps).

``` r

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
#>   Fbp omissionRate SBI_tp SBI_cr SBI_bs SBI_ps SBI_ad  SBI_m SEDI ORSS AUC
#> 1   2            0      1      1      1      1 0.9223 0.9991    0    1   1
#>      COR Spec Sens Kappa PCC TSS PRG   MAE   BIAS noPresencePoints
#> 1 0.9282    1    1     1   1   1   1 0.205 0.0317                3
```

## Presence and background as `sf` points

[`calculateMetrics()`](https://envima.github.io/poEvaluationMetrics/reference/calculateMetrics.md)
extracts values at presence and absence/background points, then calls
[`indexCalculation()`](https://envima.github.io/poEvaluationMetrics/reference/indexCalculation.md).

``` r

pres_coords <- data.frame(x = c(0.2, 0.35, 0.6), y = c(0.25, 0.7, 0.45))
bg_coords   <- data.frame(x = c(0.1, 0.5, 0.85), y = c(0.55, 0.15, 0.75))

presence <- st_as_sf(pres_coords, coords = c("x", "y"), crs = crs(r))
absence  <- st_as_sf(bg_coords, coords = c("x", "y"), crs = crs(r))

round(calculateMetrics(r, presence, absence), 4)
#> Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L = G$L,
#> : Iteration limit reached without full convergence - check carefully
#> Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L = G$L,
#> : Iteration limit reached without full convergence - check carefully
#>      Fbp omissionRate  SBI_tp SBI_cr SBI_bs SBI_ps SBI_ad  SBI_m   SEDI ORSS
#> 1 0.6667       0.6667 -0.4069 0.0657 0.0669 -0.367  0.054 0.0741 0.8774    1
#>      AUC    COR Spec   Sens  Kappa    PCC    TSS  PRG    MAE   BIAS
#> 1 0.3333 -0.173    1 0.3333 0.3333 0.6667 0.3333 -0.5 0.5495 0.1787
#>   noPresencePoints
#> 1                3
```

## Multiple evaluation designs with `performanceEstimation()`

[`performanceEstimation()`](https://envima.github.io/poEvaluationMetrics/reference/performanceEstimation.md)
combines presence–absence, presence–background (replicated background
sampling), and/or presence–artificial-absence (AOA) evaluation. It
returns one **`data.frame`** with a **`scenario`** column (`PA`, `PBG`,
`PAA`).

The chunk below keeps replication small for illustration. Turning on
`aa = TRUE` uses **CAST** and is slower; here we disable PAA and use
only PBG so the vignette builds quickly.

``` r

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
#> Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L = G$L,
#> : Iteration limit reached without full convergence - check carefully
#> Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L = G$L,
#> : Iteration limit reached without full convergence - check carefully
#> Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L = G$L,
#> : Iteration limit reached without full convergence - check carefully
#> Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L = G$L,
#> : Iteration limit reached without full convergence - check carefully

pe_result[, 1:8]
#>   scenario       Fbp omissionRate     SBI_tp     SBI_cr     SBI_bs     SBI_ps
#> 1       PA 0.6666667    0.6666667 -0.4102622 0.06091446 0.06212713 -0.3693533
#> 2      PBG 0.3095238    0.6666667 -0.1772434 0.05673406 0.05899522 -0.1458808
#>       SBI_ad
#> 1 0.04270814
#> 2 0.03601108
```

See
[`?performanceEstimation`](https://envima.github.io/poEvaluationMetrics/reference/performanceEstimation.md)
for all arguments and design options.

## Further reading

- Function reference in
  [`help(package = "poEvaluationMetrics")`](https://rdrr.io/pkg/poEvaluationMetrics/man).
- README dependency notes for **prg** and **mecofun** if you install
  from Git.
