## Overview

`dda`, Direction Dependence Analysis (DDA) provides framework for analyzing competing linear models. A target model `y ~ x` is compared to an alternate (causally reversed) model `x ~ y` through a series of diagnostic tests. DDA framework supports causal model exploration and potential confounding detection through diagnostics with higher-order moments.

* `cdda.indep()` conditional (moderation) independence property tests, including non‐linear correlation tests, Breusch–Pagan homoscedasticity tests, and the HSIC test 
* `cdda.vardist()` conditional (moderation) variable distribution‐based tests, including D'Agostino and Anscombe–Glynn tests and bootstrap CIs on higher moment differences
* `dda.indep()` independence property tests, including non‐linear correlation tests, Breusch–Pagan homoscedasticity tests, and the HSIC test 
* `dda.resdist()` residual distribution tests, including D'Agostino and Anscombe–Glynn tests and bootstrap CIs on higher moment differences
* `dda.vardist()` variable distribution‐based tests, including D'Agostino and Anscombe–Glynn tests and bootstrap CIs on higher moment differences

Conditional (moderation) functions can be plotted as well where different moderator levels are visually displayed. 

If you are new to Direction Dependence Analysis (DDA) concepts, the best place to start is the [Direction Dependence in Statistical Modeling: Methods of Analysis](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119523024) text.

## Installation

The `dda` development version can be installed from GitHub: 

```
remotes::install_github("wwiedermann/dda")
```

## Usage

```
library(dda)
```


```
n <- 1000

### generate moderator
z <- sort(rnorm(n))
z1 <- z[z <= 0]; z2 <- z[z > 0]

### x -> y when m <= 0
x1 <- rchisq(length(z1), df = 4) - 4
e1 <- rchisq(length(z1), df = 3) - 3
y1 <- 0.5 * x1 + e1

### y -> x when m > 0
y2 <- rchisq(length(z2), df = 4) - 4
e2 <- rchisq(length(z2), df = 3) - 3
x2 <- 0.25 * y2 + e2

y <- c(y1, y2); x <- c(x1, x2)
dat <- data.frame(x,y,z)

m <- lm(y ~ x*z, data = dat)
##summary(m)
```

```
mean.indep <- cdda.indep(m, pred = "x", mod = "z", data = dat, nlfun = 2,
                          modval = "mean", diff = TRUE, hetero = TRUE)

summary(mean.indep, hsic.diff = TRUE, dcor.diff = TRUE, mi.diff = TRUE)
plot.cddaindep(mean.indep, stat = "hsic.diff")
```

```
point.vardist <- cdda.vardist(m, pred = "x", mod = "z", data = dat,
                          modval = c(-1, 0, 1))

summary(point.vardist, coskew = TRUE, cokurt = TRUE)
plot(mean.vardist, stat = "rhs", ylim = c(-0.2, 0.3))
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal reproducible example on GitHub. For questions and other discussion, please contact the package maintainer.
