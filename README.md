# logsplineCDA

This folder is a lightweight R package wrapping the helper functions used in
Rozen & Canales (2025), *A Logspline-Based Approach for Analyzing Left-Censored
Concentration Data*.

## Install (local)

From an R session:

```r
# from the parent directory containing `logsplineCDA/`
install.packages("logsplineCDA", repos = NULL, type = "source")

# or, with devtools:
# devtools::install_local("logsplineCDA")
```

## Main estimation functions

All estimators take:

- `censored_sim`: numeric vector containing observations, with censored values
  recorded at their detection limit(s)
- `censored_bool`: logical vector of the same length (`TRUE` = censored)

They return a named list with `mean`, `sd`, `gm`, and `gsd`.

Example:

```r
library(logsplineCDA)

x  <- c(0.5, 0.5, 0.7, 1.2, 2.0)          # two censored at LOD=0.5
cen <- c(TRUE, TRUE, FALSE, FALSE, FALSE)

# Logspline estimators
fit <- spline_fit(x, cen)
spline_estimate(x, cen, fit)
robust_spline(x, cen, fit)

# Alternatives
ROS_estimate(x, cen)
MLE_estimate(x, cen, method = "mle")
KM_estimate(x, cen)
substitution_estimate(x, cen, method = "half")
```

## Dependencies

This package expects the following CRAN packages to be installed:

- `logspline`
- `EnvStats`
- `cubature`
- `survival`

A small subset of NADA-derived code is bundled for internal use (ROS and KM helpers).
