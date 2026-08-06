# fmrireg compatibility interface

Stable entry points to help upstream packages reuse fmriAR whitening
without rewriting existing pipelines.

## Usage

``` r
compat
```

## Format

An object of class `list` of length 6.

## Value

A list environment containing compatibility functions:

- `plan_from_phi`: Create whitening plan from AR coefficients

- `whiten_with_phi`: Apply whitening given AR coefficients

- `update_plan`: Update existing plan with new residuals

- `plan_info`: Extract information from a plan object

- `whiteness_score`: Compute whiteness metric from residuals

- `afni_restricted_plan`: Build AFNI-style restricted AR plan from root
  parameters (advanced; internal helper exposed via compat)

## Examples

``` r
# Create compatibility interface
compat_funcs <- compat

# Example: Create whitening plan from AR coefficients
phi <- c(0.3, 0.1)  # AR(2) coefficients
plan <- compat_funcs$plan_from_phi(phi, exact_first = TRUE)
#> Warning: no non-missing arguments to max; returning -Inf

# Example: Compute whiteness score
resid <- matrix(rnorm(100 * 10), 100, 10)
score <- compat_funcs$whiteness_score(resid)
```
