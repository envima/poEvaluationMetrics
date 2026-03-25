# Generate Background Points for SDM Evaluation

Generates random background points for presence-background (PBG)
evaluation of species distribution models. Points are sampled from the
environmental raster layers and returned as an `sf` object.

## Usage

``` r
generateBackgroundPoints(rasters, nPoints)
```

## Arguments

- rasters:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object containing environmental covariates.

- nPoints:

  Integer. Number of background points to generate.

## Value

An `sf` object with randomly sampled background points.
