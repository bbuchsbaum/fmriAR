# Build an AFNI-style restricted AR plan from root parameters

Build an AFNI-style restricted AR plan from root parameters

## Usage

``` r
afni_restricted_plan(
  resid,
  runs = NULL,
  parcels = NULL,
  p = 3L,
  roots,
  estimate_ma1 = TRUE,
  exact_first = TRUE
)
```

## Arguments

- resid:

  (n x v) residual matrix (used only if estimate_ma1=TRUE)

- runs:

  integer vector length n (optional)

- parcels:

  integer vector length v (optional; if provided, plan pooling='parcel')

- p:

  either 3 or 5

- roots:

  either a single list with elements named as needed - for p=3: list(a,
  r1, t1, vrt = 1.0) - for p=5: list(a, r1, t1, r2, t2, vrt = 1.0) or a
  named list of such lists keyed by parcel id (character) for per-parcel
  specs.

- estimate_ma1:

  logical, if TRUE estimate MA(1) on AR residuals to mimic AFNI's
  additive white

- exact_first:

  apply exact AR(1) scaling at segment starts (harmless here; default
  TRUE)

## Value

An `fmriAR_plan` with `method = "afni"` that can be supplied to
[`whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.md).

## Examples

``` r
NULL
#> NULL
```
