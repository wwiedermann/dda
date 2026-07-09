# dda 1.1.1

---

### New features

- Native C++ (Rcpp) backend for the Hilbert-Schmidt Independence Criterion in
  `src/hsic.cpp`, replacing the previous pure-R implementation. Exposes
  `hsic()`, `hsic_test()` (methods `"gamma"`, `"permutation"`, `"eigenvalue"`,
  `"bootstrap"`), and `hsic_resid_test()`.
- `dda.bagging()`: bootstrap-aggregated DDA to assess the stability and
  robustness of direction-of-dependence decisions, with `print`/`summary`
  methods.
- `robust` argument for `dda.indep()` and `dda.resdist()`: Siegel (1982)
  repeated-median estimation for the competing models, with a companion robust
  Breusch-Pagan test (`bptestrobust`).

### Internal

- HSIC C++ routines migrated from an inline `Rcpp::sourceCpp()` call to a
  compiled `src/` backend registered via `useDynLib`.
- Removed the superseded `boot_hsic_test.R`; consolidated the duplicated
  `nlcor.test` into a single definition.


# dda 0.1.1

## dda 0.1.1

---

### Bug fixes

-	In `dda.indep`, the bootstrap HSIC method argument is now `hsic.method = "bootstrap"` (previously `hsic.method = "boot"` in 0.1.0).
-	`dda.vardist` and `dda.resdist` now consistently return error messages when the number of bootstrap replications (B) is too small for `boot.type = "bca"` (previously, behavior could be inconsistent).
-	Various minor documentation changes and clarifications.


## dda 0.1.0

---

### Initial release

- First CRAN release of `dda`.
- Includes five core `dda` functions and S3 generics where applicable (`print`, `summary`, `plot`).
- Documentation provided for all exported/user-facing functions.

## version 0.0.0.9000

---

### Development

- Added NEWS.md using [newsmd](https://github.com/Dschaykib/newsmd) package.
