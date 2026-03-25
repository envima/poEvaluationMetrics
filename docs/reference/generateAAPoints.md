# Generate Artificial Absence Points for SDM Evaluation

Generates artificial absence (AA) points for presence-artificial-absence
(PAA) evaluation of species distribution models. Uses the Area of
Applicability (AOA) mask to sample points from regions environmentally
dissimilar to presence locations.

## Usage

``` r
generateAAPoints(aa_mask, nPoints)
```

## Arguments

- aa_mask:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object containing environmentally not suitable areas.

- nPoints:

  Integer. Number of artificial absence points to generate.

## Value

An `sf` object with sampled artificial absence points.
