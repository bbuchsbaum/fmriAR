# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

fmriAR is an R package that provides fast AR/ARMA prewhitening for fMRI design and data matrices. It uses C++ (via Rcpp and RcppArmadillo) for performance-critical computations.

## Development Commands

### Build and Load
```r
# Compile C++ attributes
Rcpp::compileAttributes()

# Generate documentation
devtools::document()

# Load package for development
devtools::load_all()
```

### Testing
```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-whitening.R")
```

### Package Check
```r
# Full R CMD check
devtools::check()
```

## Architecture

### Core Components

1. **Noise Model Estimation** (`R/fit_and_apply.R`)
   - `fit_noise()`: Estimates AR/ARMA parameters from residuals
   - Supports run-aware and censor-aware fitting
   - Auto-selects model order when `p = "auto"`
   - Pooling options: "global", "run", or "none"

2. **Whitening Application** (`R/fit_and_apply.R`)
   - `whiten_apply()`: Applies prewhitening transformation to design matrix X and data Y
   - Returns whitened versions suitable for standard linear model fitting
   - Handles multi-run fMRI data with proper boundary conditions

3. **C++ Backend** (`src/`)
   - `arma_whiten_inplace()`: Core whitening implementation using RcppArmadillo
   - OpenMP support (optional) for parallel processing
   - In-place operations for memory efficiency

4. **Statistical Utilities**
   - `sandwich_from_whitened_resid()` (`R/sandwich.R`): Computes robust standard errors
   - `acorr_diagnostics()` (`R/acorr.R`): Autocorrelation diagnostics for residuals
   - `pacf_cpp()` (`R/pacf_helpers.R`): Partial autocorrelation computations

### Typical Workflow

The package is designed for a specific fMRI preprocessing workflow:

1. Fit initial OLS model to get residuals
2. Estimate noise model from residuals using `fit_noise()`
3. Apply whitening to both design matrix and data using `whiten_apply()`
4. Fit final model on whitened data
5. Compute robust standard errors and diagnostics

## Build Requirements

- R >= 4.0
- C++14 compiler
- RcppArmadillo (linked automatically)
- Optional: OpenMP support (uncomment in `src/Makevars` for parallel processing)