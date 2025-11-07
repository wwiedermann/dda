## dda 0.1.1

---


### Error Patch

- In `dda.indep` for the bootstrap HSIC method, now specify `hsic.method = "bootstrap"`.
- `dda.vardist` and `dda.resdist` should now consistently return errors when bootstraps (`B`) are too low (specifically for `boot.type = "bca"`)



## dda 0.1.0

---


### Initial Release

- First CRAN release of `dda`.
- Includes five core `dda` functions along with `print`, `summary`, and `plot` generics when applicable.
- Documentation provided for all external functions.


## version 0.0.0.9000

---

### NEWS.md setup

- added NEWS.md creation with [newsmd](https://github.com/Dschaykib/newsmd)

