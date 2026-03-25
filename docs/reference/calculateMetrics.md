# Calculate Performance Metrics for a single replicate/model/map

Calculates a suite of performance metrics for a single replicate of SDM
evaluation. Supports presence-background (PBG),
presence-artificial-absence (PAA), and presence-absence (PA) datasets.

## Usage

``` r
calculateMetrics(prediction, presence, absence_or_bg_sf)
```

## Arguments

- prediction:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object containing model predictions.

- presence:

  An `sf` object of known presence locations.

- absence_or_bg_sf:

  An `sf` object of absence, background, or artificial absence points.

## Value

A one-row `data.frame` with calculated metrics for the replicate.
