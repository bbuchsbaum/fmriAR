# Pretty-print an fmriAR whitening plan

Pretty-print an fmriAR whitening plan

## Usage

``` r
# S3 method for class 'fmriAR_plan'
print(x, ...)
```

## Arguments

- x:

  An object returned by
  [`fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.md).

- ...:

  Unused; included for S3 compatibility.

## Value

The input plan, invisibly.

## Examples

``` r
resid <- matrix(rnorm(60), 20, 3)
plan <- fit_noise(resid, method = "ar", p = 2)
print(plan)
#> fmriAR whitening plan
#>   Method: AR
#>   Orders: p = 2, q = 0
#>   Pooling: global
#>   Exact first-sample scaling: AR(1)
#>   Coefficients:
#>     global: phi = 0.180, -0.133
```
