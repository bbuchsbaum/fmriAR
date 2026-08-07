# CRAN submission: fmriAR 0.3.3

Update to the CRAN version 0.3.2. User-facing changes are listed in NEWS.md.

## R CMD check results

0 errors | 0 warnings | 0 notes

Local `R CMD check --as-cran --no-manual` on macOS (aarch64-apple-darwin20,
R 4.5.1) additionally shows one warning that does not originate in package
code: an external compiler message from `R_ext/Boolean.h` in the R headers on
this local macOS/Homebrew clang toolchain (`unknown warning group
'-Wfixed-enum-extension'`). It was also present, and accepted, in the 0.3.2
submission.

Also checked on win-builder (R-devel, R-release), the macOS builder, and
GitHub Actions (Linux, macOS, Windows).

## Reverse dependencies

There are no reverse dependencies on CRAN.
