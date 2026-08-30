# dda 0.2.0

---

### New features

- `dda.bagging()`: bootstrap-aggregated DDA to assess the stability and
  robustness of direction-of-dependence decisions across resamples, with
  `print`/`summary` methods and aggregated OLS summaries
  (`print_ols_summary()`, `summary_ols()`). p-values are combined using the
  harmonic mean p-value approach (Wilson, 2019). Aggregation statistics
  include `"mean"`, `"median"`, `"trimmed"`, `"winsorized"`, `"midhinge"`,
  and `"tukey"`; `inner_B` caps the inner resampling budget of each
  per-iteration DDA call.
- `robust` argument for `dda.indep()`, `cdda.indep()`, and `dda.resdist()`:
  Siegel (1982) repeated-median estimation for the causally competing
  models, with a companion robust Breusch-Pagan test (`bptestrobust`).

### Bug fixes

-	In `dda.indep`, `boot.type = "bca"` no longer aborts when the
	acceleration constant cannot be calculated; the function now warns and
	falls back to percentile intervals.
-	Various minor documentation changes and clarifications.

### Internal

- Distance correlation difference statistics are now computed with
  `dccpp::dcor()`. HSIC inference continues to be supplied by
  `dHSIC::dhsic.test()`.
- Removed the superseded `boot_hsic_test.R`; consolidated the duplicated
  `nlcor.test` into a single definition.
- Examples with non-trivial run time now pair a fast, low-`B` illustration
  with a `\dontrun{}` block showing realistic resampling settings.


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
