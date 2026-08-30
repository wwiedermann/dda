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
* In dda.indep(), boot.type = "bca" now warns and falls back to percentile
  intervals when the acceleration constant cannot be calculated, rather than
  aborting.
* Examples with non-trivial run time pair a fast, low-B illustration with a
  \dontrun{} block showing realistic resampling settings.
* See NEWS.md for the full changelog.

## Note on dependencies

This release remains pure R (NeedsCompilation: no). HSIC inference is provided
by dHSIC, which has been added back to Imports; distance correlation difference
statistics use dccpp. There is no compiled code in this package.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.5.2 -- 0 errors | 0 warnings | 0 notes
* TODO before submitting -- confirm on:
    - win-builder (R-devel and R-release): devtools::check_win_devel() / check_win_release()
    - macOS (R-release): devtools::check_mac_release()
    - R-hub (Linux; gcc + clang): rhub::rhub_check()
* TODO before submitting -- restore the maintainer and author e-mail addresses
  in DESCRIPTION (and man/dda-package.Rd) to
  Wolfgang Wiedermann <wiedermannw@missouri.edu> and
  Megan Hirni <mjhirni@outlook.com>. They are temporarily set to
  mj.hirni@missouri.edu so that check notifications route to the second author
  while pre-submission checks are running.

## Reverse dependencies

There are no reverse dependencies on CRAN. Confirmed with
revdepcheck::revdep_check() on 2026-08-27: 0 reverse dependencies
checked, 0 new problems.

## Maintainer

Wolfgang Wiedermann <mj.hirni@missouri.edu>
