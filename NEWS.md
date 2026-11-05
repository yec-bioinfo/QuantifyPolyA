# QuantifyPolyA 0.4.0

## Breaking / behaviour changes

* `Remove.IP()` and `is.internal.priming()` gain three arguments controlling
 their internal priming detection: `flank_len` (default `15`), `win_size`
 (default `10`) and `min_A` (default `8`). The detection window is now
 `coord +/- flank_len` and the threshold is `min_A` adenines in a sliding
 window of `win_size` bases, replacing the previously hard-coded
 `coord +/- 10` / `>=7` / window `10`.
* Because the defaults changed, results from `Remove.IP()` are not identical
 to version 0.3.0. To reproduce the old behaviour use
 `Remove.IP(QpolyA, fasta, flank_len = 10, win_size = 10, min_A = 7)`.

## Dependency changes

* `makeTxDbFromGFF()` moved from GenomicFeatures to the **txdbmaker** package
 in Bioconductor 3.19 and was removed from GenomicFeatures in 3.23.
 `polyAsite.annotation()` now calls `txdbmaker::makeTxDbFromGFF()`, and
 `txdbmaker` is a new dependency.
* All packages moved from `Depends` to `Imports`. `library(QuantifyPolyA)` no
 longer attaches DESeq2, ggplot2, dplyr and the other dependencies to the
 search path.
* `rtracklayer` is no longer a dependency; it is only needed to import a GTF
 for the examples in the vignette and is now listed under `Suggests`.
* `utils`, `stats`, `methods` and `grDevices` are now imported explicitly.

## Package infrastructure

* Added `Authors@R`, `URL`, `BugReports`, `Remotes` and `VignetteBuilder`.
* Added a package vignette based on the SourceForge user manual.
* Added `README.md`, `NEWS.md` and GitHub Actions workflows for
 `R CMD check` and pkgdown.
* Removed the top-level `if (!require(...)) install.packages(...)` calls from
 the package source; installing packages during `R CMD INSTALL` is not
 allowed and these calls never ran for users installing the built package.
* Version numbering: the package now follows the usual
 major.minor.patch scheme (`0.4.0`).

# QuantifyPolyA 0.3.0

* Last version distributed through SourceForge.
