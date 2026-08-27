## Submission - dda 0.2.0

This is a feature update to the CRAN package 'dda' (current CRAN version 0.1.1).
All changes are backward-compatible.

## Summary of changes

* New C++ (Rcpp) backend for the Hilbert-Schmidt Independence Criterion
  (src/hsic.cpp), replacing the previous pure-R implementation. New exported
  functions hsic() and hsic.test(), with inference methods "gamma",
  "permutation", "eigenvalue", and "bootstrap".
* New dda.bagging(): bootstrap-aggregated Direction Dependence Analysis for
  assessing the stability and robustness of direction-of-dependence decisions,
  with accompanying print() and summary() methods.
* New 'robust' argument for dda.indep() and dda.resdist(): Siegel (1982)
  repeated-median regression, with a companion robust Breusch-Pagan test.
* See NEWS.md for the full changelog.

## Note on compiled code

Unlike the previous release (0.1.1, pure R), this version contains compiled
C++ code in src/, interfaced via Rcpp (NeedsCompilation: yes). This is the main
change reviewers may wish to note.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.5.2 -- 0 errors | 0 warnings | 0 notes
* TODO before submitting -- confirm on:
    - win-builder (R-devel and R-release): devtools::check_win_devel() / check_win_release()
    - macOS (R-release): devtools::check_mac_release()
    - R-hub (Linux; gcc + clang, and a sanitizer flavour for the new C++):
      rhub::rhub_check()

## Reverse dependencies

There are no reverse dependencies on CRAN. Confirmed with
revdepcheck::revdep_check() on 2026-08-27: 0 reverse dependencies
checked, 0 new problems.

## Maintainer

Wolfgang Wiedermann <wiedermannw@missouri.edu>
