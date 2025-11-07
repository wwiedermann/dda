## dda 0.1.1

---

### Bug fixes

- In `dda.indep` the bootstrap HSIC method argument is now `hsic.method = "bootstrap"` (previously `hsic.method = "boot"` in 0.1.0).
- `dda.vardist` and `dda.resdist` now consistently throw an error when the number of bootstrap replications (`B`) is too low for `boot.type = "bca"` (previously behavior could be inconsistent).

## dda 0.1.0

---

### Initial release

- First CRAN release of `dda`.
- Includes five core `dda` functions and S3 generics where applicable (`print`, `summary`, `plot`).
- Documentation provided for all exported/user-facing functions.

## version 0.0.0.9000

---

### Development

- Added automatic NEWS.md creation using [newsmd](https://github.com/Dschaykib/newsmd).
