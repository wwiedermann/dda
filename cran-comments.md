# cran-comments.md for dda (0.1.1)

### R CMD check results

0 errors | 0 warnings | 0 notes

We saw 3 new problems:

* dda.indep
* dda.vardist
* dda.resdist

Both maintainers were notified on Oct 29 (~2 weeks ago) and supplied with patches.

Package: dda
Version: 0.1.1
Date: 2025-11-07

Summary of changes
- Bug fixes: Correction in the  HSIC method argument in `dda.indep` (now `hsic.method = "bootstrap"`), and `dda.vardist` and `dda.resdist` now consistently return error messages when the number of bootstrap replications (B) is too small for `boot.type = "bca"`. 
- See NEWS.md for the changelog.

Check results
- Ran R CMD check –as-CRAN locally (Windows 11 and macOS Ventura) and on win-builder. - Tested on the current CRAN release platform(s) and on Windows via the win-builder service.
- Tested on the current CRAN release platform(s) and on Windows via the win-builder service. 

Package details for CRAN reviewers
- The package contains only R code (no compiled code or use of src/).  
- No special SystemRequirements (no external software required).
- Examples, tests, and vignettes should not require internet access or non-standard hardware.
- All user-facing (exported) functions are documented; internal helper functions are marked internal.

Contact
- Maintainer: Wolfgang Wiedermann (wiedermannw@missouri.edu)
- If you need additional details, test logs, or versions of R used for testing, please reach out via email.
