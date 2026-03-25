# Calculate SDM Metrics

Computes a comprehensive set of species distribution model (SDM)
evaluation metrics including correlation, AUC, PRG, MAE, BIAS, and
several skill scores.

## Usage

``` r
indexCalculation(inputDF, prediction)
```

## Arguments

- inputDF:

  A data frame with columns `predicted` and `observed`.

- prediction:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object containing predicted values.

## Value

A one-row `data.frame` containing SDM evaluation metrics (see Details).

## Details

Columns include correlation `COR`; `omissionRate`; Boyce-related
Spearman indices `SBI_tp`, `SBI_cr`, `SBI_bs`, `SBI_ps`, `SBI_ad`,
`SBI_m`; `Fbp`, `SEDI`, `ORSS`; `AUC`, `PRG`, `MAE`, `BIAS`;
[`mecofun::evalSDM`](https://rdrr.io/pkg/mecofun/man/evalSDM.html)
outputs `TSS`, `Kappa`, `PCC`, `Sens`, `Spec`; and `noPresencePoints`.
