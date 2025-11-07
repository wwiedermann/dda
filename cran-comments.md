# cran-comments.md for dda (0.1.1)

Package: dda
Version: 0.1.1
Date: 2025-11-07

Summary of changes
- Bug fixes only: corrected the bootstrap HSIC method name in dda.indep (now `hsic.method = "bootstrap"`), and made `dda.vardist` and `dda.resdist` consistently error when the number of bootstrap replications (`B`) is too low for `boot.type = "bca"`.
- See NEWS.md for the changelog.

Check results
- Ran R CMD check --as-cran locally (Windows 11 and macOS Sequoia) plus on win-builder.
- Tested on the current CRAN release platform(s) and on Windows via the win-builder service. 

Package details for CRAN reviewers
- The package contains only R code (no compiled code or use of src/).  
- No special SystemRequirements (no external software required).
- Examples, tests, and vignettes should not require internet access or non-standard hardware.
- All user-facing (exported) functions are documented; internal helper functions are marked internal.

Contact
- Maintainer: Wolfgang Wiedermann (wiedermannw@missouri.edu)
- If you need additional details, test logs, or versions of R used for testing, please let me know and I will provide them.
