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
  specs. As in AFNI, `a`, `r1`, `r2` must be non-negative (values above
  0.95 are clamped to 0.95) and angles are clamped to `[0, pi]`. `vrt`
  is AFNI's signal-to-total variance ratio `s^2 / (s^2 + w^2)` for
  additive white noise `w`; it is used when `estimate_ma1 = FALSE`.

- estimate_ma1:

  logical. If `TRUE`, estimate an MA(1) term from the AR-filtered voxels
  to stand in for AFNI's additive white noise (`vrt` is then ignored).
  If `FALSE`, the white-noise share is taken from `vrt`: for `vrt < 1`
  the plan is the exact ARMA(p, p) equivalent of AR(p) plus white noise,
  so its autocorrelation is `vrt` times the AR autocorrelation at every
  non-zero lag; `vrt = 1` (default) gives the pure AR(p).

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
