# Symmetric Extremal Dependence Index (SEDI)

Evaluates prediction skill for rare events.

## Usage

``` r
sedi(actual, predicted, eps = 1e-10)
```

## Arguments

- actual:

  Vector of observed binary outcomes.

- predicted:

  Vector of predicted binary outcomes.

- eps:

  Small number to prevent division by zero. Default = 1e-10.

## Value

Numeric SEDI value.
