# Top-q% Capture Rate

Computes proportion of presences captured in top q% of predicted values.

## Usage

``` r
tcr(predicted, observed, prediction, q = 0.1)
```

## Arguments

- predicted:

  Vector of predicted values.

- observed:

  Vector of observed binary outcomes.

- prediction:

  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object used to compute threshold.

- q:

  Numeric. Top proportion threshold (default 0.10 = top 10%).

## Value

Numeric proportion of presences captured.
