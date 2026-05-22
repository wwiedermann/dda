# cran-comments.md for dda (0.1.2)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Issues addressed

Two main changes occur in this release:

- The primary author (Wiedermann) was notified on 30 April 2026
of the potential `dHSIC` CRAN package archival. This patch release provides
a local `hsic` function to replace the dependency on `dHSIC`.

- **dda.indep** now utilizes the local `hsic` function which computes the
Hilbert-Schmidt Independence Criterion (HSIC) for two variables. 

- **hsic** contains options for different kernel types and 
methods for estimating the null distribution.


## Test environments

- Local: Windows 11, R release
- Local: macOS Ventura, R release
- win-builder: R release and R devel

## Package details

- Contains only R code; no compiled code and no `src/` directory.
- No external `SystemRequirements`.
- All examples, tests, and vignettes run without internet access or
  non-standard hardware.
- All exported functions are documented; internal helpers are marked
  with `@keywords internal`.

## Changes in this version

See `NEWS.md` for the changelog.

## Contact

Maintainer: Wolfgang Wiedermann <wiedermannw@missouri.edu>
Additional test logs or platform details are available on request.
