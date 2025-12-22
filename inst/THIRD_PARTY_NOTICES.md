# Third-party notices

This package bundles or derives code from the following third-party R packages.

## NADA

Portions of the files in `R/00_nada_All.R`, `R/01_nada_cen.R`, `R/02_nada_km.R`, and
`R/03_nada_ros.R` are copied from the CRAN package **NADA**.

- CRAN page: https://CRAN.R-project.org/package=NADA
- License: GPL (>= 2) (listed on CRAN as GPL-2 | GPL-3)

## logspline

The function `oldlogspline()` in `R/10_lod_methods.R` is a lightly edited copy
of `logspline::oldlogspline()` (Kooperberg and Stone, 1992 algorithm). The
function still calls the compiled routine shipped with **logspline**.

- CRAN page: https://CRAN.R-project.org/package=logspline
- License: Apache License 2.0

## EnvStats

This package depends on **EnvStats** for lognormal MLE and ROS-style imputation.

- CRAN page: https://CRAN.R-project.org/package=EnvStats
- License: GPL (>= 3)

## cubature

This package depends on **cubature** for numerical integration.

- CRAN page: https://CRAN.R-project.org/package=cubature
- License: GPL-3
