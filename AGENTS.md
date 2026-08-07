# AGENTS.md

`fmriAR` is an R package (with an Rcpp / RcppArmadillo C++ backend in `src/`) for
fast AR/ARMA prewhitening of fMRI design and data matrices. It has no GUI or
long-running service — the "application" is the R library itself, exercised from
an R session or `Rscript`.

Standard developer commands (build, document, test, check) are documented in
`CLAUDE.md`. Prefer those; the notes below only add cloud/environment-specific
context that is not obvious from the existing docs.

## Cursor Cloud specific instructions

- R (4.6.x) is installed system-wide via apt from the CRAN repo, and R package
  dependencies live in the per-user library at
  `~/R/x86_64-pc-linux-gnu-library/4.6` (R picks this up automatically; no env
  var needed). The startup update script only refreshes R package deps, so R
  itself, `gfortran`, and `pandoc` come from the VM snapshot.
- Package installs use the Posit Package Manager Linux binary repo
  (`https://packagemanager.posit.co/cran/__linux__/noble/latest`) for speed. A
  matching `HTTPUserAgent` must be set for binaries to be served; installing from
  source (plain `install.packages` with the default CRAN mirror) still works but
  compiles everything and is much slower.
- Loading the package compiles the C++ backend: run
  `Rscript -e 'devtools::load_all(".")'`. After changing any `src/*.cpp` file,
  regenerate bindings with `Rscript -e 'Rcpp::compileAttributes()'` before
  reloading.
- `devtools::check()` / `devtools::document()` here run roxygen2 8.1.0, which is
  newer than the `RoxygenNote: 7.3.3` the repo was built with. Documenting will
  rewrite `man/*.Rd` and bump `RoxygenNote` in `DESCRIPTION`. Do NOT commit those
  regenerated doc/DESCRIPTION churn unless you intentionally mean to update
  roxygen output — revert with `git checkout -- DESCRIPTION man/` after a check.
- `devtools::check()` reports exactly one benign `NOTE`: "Compilation used the
  following non-portable flag(s): `-mno-omit-leaf-frame-pointer`". This flag comes
  from R's own system `CFLAGS` on this Ubuntu build, not from the package, and can
  be ignored (0 errors / 0 warnings otherwise).
- OpenMP is intentionally disabled by default (commented out in `src/Makevars`),
  so the `#pragma omp ... [-Wunknown-pragmas]` compile warning from
  `fmriAR_whiten.cpp` is expected and harmless.
