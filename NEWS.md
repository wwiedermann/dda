# dda 0.1.2

## dda 0.1.2

---

### HSIC function added

- Added `hsic` function to compute the Hilbert-Schmidt Independence Criterion (HSIC) for two variables, with options for different kernel types and methods for estimating the null distribution.

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
