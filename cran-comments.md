## Submission - dda 0.2.0

This is a feature update to the CRAN package 'dda' (current CRAN version 0.1.1).
All changes are backward-compatible.

## Summary of changes

* New dda.bagging(): bootstrap-aggregated Direction Dependence Analysis for
  assessing the stability and robustness of direction-of-dependence decisions,
  with accompanying print() and summary() methods and aggregated OLS summaries.
  p-values are combined across resamples using the harmonic mean p-value.
* New 'robust' argument for dda.indep(), cdda.indep(), and dda.resdist():
  Siegel (1982) repeated-median regression, with a companion robust
  Breusch-Pagan test.
* See NEWS.md for the full changelog.

## Note on dependencies

This release remains pure R (NeedsCompilation: no). Distance correlation difference
statistics use dccpp. There is no compiled code in this package.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.5.2
* win-builder: R-devel and R-release
* macOS: R-release

All test environments returned 0 errors | 0 warnings | 0 notes.

## Reverse dependencies

There are no reverse dependencies on CRAN. Confirmed with
revdepcheck::revdep_check() on 2026-09-02: 0 reverse dependencies
checked, 0 new problems.

## Maintainer

Wolfgang Wiedermann <wiedermannw@missouri.edu>
