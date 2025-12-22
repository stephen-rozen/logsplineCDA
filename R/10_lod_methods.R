# logsplineCDA: user-facing estimators for left-censored (LOD) concentration data
#
# This file contains the core analysis functions intended for end users.
# Simulation-only helpers are intentionally excluded from the package API.


# ---- Modified oldlogspline() (from logspline) ------------------------------
#
# This is a lightly edited copy of logspline::oldlogspline(). The main edits
# are:
# * noisy console output moved from cat/print into message() where feasible.
# * allows use of suppressMessages() in simulation settings.
#
# NOTE: This function still calls the compiled routine in logspline.
#
# (Original algorithm: Kooperberg and Stone, 1992)

oldlogspline <- function(uncensored, right, left, interval, lbound, ubound,
                         nknots, knots, penalty, delete = TRUE) {
  nsample <- rep(0, 6)
  if (!missing(uncensored)) uncensored <- logspline::unstrip(uncensored)
  if (!missing(right))      right      <- logspline::unstrip(right)
  if (!missing(left))       left       <- logspline::unstrip(left)
  if (!missing(interval))   interval   <- logspline::unstrip(interval)
  if (!missing(knots))      knots      <- logspline::unstrip(knots)

  if (!missing(interval)) {
    if (length(interval[1, ]) != 2) stop("interval must have two columns")
    if (min(abs(interval[, 1] - interval[, 2])) < 0)
      stop("not all lower bounds smaller than upper bounds")
    nsample[3] <- length(interval) / 2
    nsample[1] <- length(interval) / 2
    if (!missing(lbound)) interval[interval[, 1] < lbound, 1] <- lbound
    if (!missing(ubound)) interval[interval[, 2] > ubound, 2] <- ubound
    sample <- as.vector(t(interval))
    ror <- order(interval[, 1], interval[, 2])
    if (nsample[3] > 1) {
      ro1 <- interval[ror[(1:(nsample[3] - 1))], 1] == interval[ror[2:nsample[3]], 1]
      ro2 <- interval[ror[(1:(nsample[3] - 1))], 2] == interval[ror[2:nsample[3]], 2]
      nsample[6] <- nsample[3] - sum(ro1 + ro2 == 2)
    } else nsample[6] <- 1
  }

  if (!missing(uncensored)) {
    uncensored2 <- uncensored[!is.na(uncensored)]
    u2 <- length(uncensored) - length(uncensored2)
    if (u2 > 0) message(paste("***", u2, " NAs ignored in uncensored"))
    uncensored <- uncensored2
    if (nsample[1] > 0) sample <- c(uncensored, sample)
    if (nsample[1] == 0) sample <- uncensored
    nsample[1] <- length(uncensored) + nsample[1]
    nsample[2] <- length(uncensored)
    uncensored <- sort(uncensored)
    if (nsample[2] > 1)
      nsample[6] <- sum(uncensored[2:nsample[2]] != uncensored[1:(nsample[2] - 1)]) + 1 + nsample[6]
    else nsample[6] <- nsample[6] + 1
  }

  if (nsample[1] == 0)
    stop("you either need uncensored or interval censored data")

  if (!missing(right)) {
    if (nsample[1] > 0) sample <- c(sample, right)
    if (nsample[1] == 0) sample <- right
    nsample[1] <- length(right) + nsample[1]
    nsample[4] <- length(right)
    right <- sort(right)
    if (nsample[4] > 1) {
      nsample[6] <- sum(right[2:nsample[4]] != right[1:(nsample[4] - 1)]) + 1 + nsample[6]
    } else nsample[6] <- nsample[6] + 1
  }

  if (!missing(left)) {
    if (nsample[1] > 0) sample <- c(sample, left)
    if (nsample[1] == 0) sample <- left
    nsample[1] <- length(left) + nsample[1]
    nsample[5] <- length(left)
    left <- sort(left)
    if (nsample[5] > 1) {
      nsample[6] <- sum(left[2:nsample[5]] != left[1:(nsample[5] - 1)]) + 1 + nsample[6]
    } else nsample[6] <- nsample[6] + 1
  }

  if (missing(penalty)) penalty <- log(nsample[1])
  n1 <- 4 * nsample[1]^0.2 + 1
  if (!missing(nknots)) n1 <- nknots + 1
  if (!missing(knots))  n1 <- length(knots) + 1

  if (!missing(knots)) {
    nknots <- length(knots)
    knots <- sort(knots)
    iautoknot <- 0
    if (knots[1] > min(sample)) stop("first knot must be smaller than smallest sample")
    if (knots[nknots] < max(sample)) stop("last knot should be larger than largest sample")
  } else {
    if (missing(nknots)) nknots <- 0
    knots <- vector(mode = "double", length = max(nknots, 50))
    iautoknot <- 1
  }

  xbound <- c(1, 0, 0, 0, 0)
  if (!missing(lbound)) {
    xbound[2] <- 1
    xbound[3] <- lbound
    if (lbound > min(sample)) stop("lbound should be smaller than smallest sample")
  }
  if (!missing(ubound)) {
    xbound[4] <- 1
    xbound[5] <- ubound
    if (ubound < max(sample)) stop("ubound should be larger than largest sample")
  }

  SorC <- vector(mode = "integer", length = 35)
  SorC[1] <- 1
  SorC[17] <- 0
  nsample[6] <- nsample[6] - 1

  if (length(table(sample)) < 3) stop("Not enough unique values")

  z <- .C("logcensor",
          as.integer(c(delete, 0, 0, 0, 0)),
          as.integer(c(iautoknot, 0, 0, 0, 0)),
          as.double(c(sample, 0, 0, 0, 0)),
          as.integer(c(nsample, 0, 0, 0, 0)),
          bd = as.double(c(xbound, 0, 0, 0, 0)),
          SorC = as.integer(c(SorC, 0, 0, 0, 0)),
          nk = as.integer(nknots),
          kt = as.double(c(knots, 0, 0, 0, 0)),
          cf = as.double(c(knots, 0, 0, 0, 0)),
          as.double(c(penalty, 0, 0, 0, 0)),
          as.double(c(sample, 0, 0, 0, 0)),
          as.double(c(sample, 0, 0, 0, 0)),
          logl = as.double(rep(0, n1 + 1 + 10)),
          PACKAGE = "logspline")

  SorC <- z$SorC

  if (SorC[1] == -1 && SorC[28] == 0 && nsample[1] != nsample[2] && nsample[2] > 15) {
    SorC <- vector(mode = "integer", length = 35)
    SorC[1] <- 1
    SorC[17] <- 1
    z <- .C("logcensor",
            as.integer(c(delete, 0, 0, 0, 0)),
            as.integer(c(iautoknot, 0, 0, 0, 0)),
            as.double(c(sample, 0, 0, 0, 0)),
            as.integer(c(nsample, 0, 0, 0, 0)),
            bd = as.double(c(xbound, 0, 0, 0, 0)),
            SorC = as.integer(c(SorC, 0, 0, 0, 0)),
            nk = as.integer(nknots),
            kt = as.double(c(knots, 0, 0, 0, 0)),
            cf = as.double(c(knots, 0, 0, 0, 0)),
            as.double(c(penalty, 0, 0, 0, 0)),
            as.double(c(sample, 0, 0, 0, 0)),
            as.double(c(sample, 0, 0, 0, 0)),
            logl = as.double(rep(0, n1 + 1 + 10)),
            PACKAGE = "logspline")
  }

  bound <- c(z$bd[2], z$bd[3], z$bd[4], z$bd[5])
  SorC <- z$SorC

  if (abs(SorC[1]) > 2) {
    for (i in 3:abs(SorC[1])) message(paste("===> warning: knot ", SorC[i - 1], " removed - double knot"))
    if (SorC[1] < 0) SorC[1] <- -1
    if (SorC[1] == 23) SorC[1] <- -3
  }

  if (abs(SorC[1]) > 3) {
    message("* several double knots suggests that your data is strongly rounded")
    SorC[1] <- 1
  }

  if (SorC[1] == -3) stop("* too many double knots")
  if (SorC[1] == -1 && SorC[28] == 0) stop("* no convergence")
  if (SorC[28] > 0)
    message(paste("* convergence problems, smallest number of knots tried is ", SorC[28] + 1, " *"))
  if (SorC[1] == 2)  stop("* sample is too small")
  if (SorC[1] == -2) stop(paste("* too many knots, at most ", SorC[2], "knots possible"))

  if (delete && SorC[28] > 0) delete <- 3

  coef <- z$cf[1:(z$nk + 2)]
  uu <- 3:z$nk
  if (delete == FALSE) uu <- 1

  fit <- list(
    coef    = coef,
    knots   = z$kt[1:z$nk],
    bound   = bound,
    logl    = z$logl[uu],
    penalty = penalty,
    sample  = nsample[1],
    delete  = delete
  )

  class(fit) <- "oldlogspline"
  fit
}


# ---- LOD-handling methods --------------------------------------------------

# internal helper
estimate <- function(object) {
  list(
    mean = mean(object),
    sd   = stats::sd(object),
    gm   = exp(mean(log(object))),
    gsd  = exp(stats::sd(log(object)))
  )
}


#' NADA-style ROS estimate
#'
#' Uses `ros()` (bundled, NADA-derived) and returns the summary stats of the
#' modeled/imputed values.
#'
#' @param censored_sim Numeric vector with observed values (LOD at censored indices).
#' @param censored_bool Logical vector (TRUE = censored).
#' @return Named list with mean, sd, gm, gsd.
#' @export
ROS_estimate <- function(censored_sim, censored_bool) {
  myros <- ros(obs = censored_sim, censored = censored_bool)
  ros_fit <- as.data.frame(myros)$modeled
  estimate(ros_fit)
}


#' Lognormal MLE (and related variants) for left-censored data
#'
#' Implements three estimators based on `EnvStats::elnormCensored(..., "mle")`:
#' * `"mle"`: standard lognormal MLE, reported as arithmetic mean/SD and geometric mean/SD.
#' * `"bcmle"`: bias-corrected mean based on a truncated series expansion (K=4).
#' * `"rmle"`: robust MLE imputation using Helsel-Cohn plotting positions.
#'
#' @param censored_sim Numeric vector with observed values (LOD at censored indices).
#' @param censored_bool Logical vector (TRUE = censored).
#' @param method One of `"mle"`, `"bcmle"`, `"rmle"`.
#' @return Named list with mean, sd, gm, gsd.
#' @export
MLE_estimate <- function(censored_sim, censored_bool, method = c("mle", "bcmle", "rmle")) {
  method <- match.arg(method)

  fit <- EnvStats::elnormCensored(censored_sim, censored_bool, "mle")
  meanlog <- unname(fit$parameters["meanlog"])
  sdlog   <- unname(fit$parameters["sdlog"])

  if (method == "mle") {
    return(list(
      mean = exp(meanlog + sdlog^2 / 2),
      sd   = sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)),
      gm   = exp(meanlog),
      gsd  = exp(sdlog)
    ))
  }

  if (method == "bcmle") {
    n <- length(censored_sim)
    g <- sdlog^2 / 2

    # ψ(g) with K = 4 (so total terms = 1 + 4 = 5)
    psi5 <- {
      res <- 1
      for (k in 1:4) {
        num <- (n - 1)^(2 * k - 1) * g^k
        den_prod <- if (k > 1) prod(n + seq(1, by = 2, length.out = k - 1)) else 1
        den <- n^k * den_prod * factorial(k)
        res <- res + num / den
      }
      res
    }

    mean_bc <- exp(meanlog) * psi5
    sd_bc   <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))

    return(list(
      mean = mean_bc,
      sd   = sd_bc,
      gm   = exp(meanlog),
      gsd  = exp(sdlog)
    ))
  }

  # robust MLE: impute censored values at fitted quantiles
  pp <- hc.ppoints(censored_sim, censored_bool)
  imputed_values <- stats::qlnorm(pp[censored_bool], meanlog, sdlog)
  all_vals <- censored_sim
  all_vals[censored_bool] <- imputed_values
  estimate(all_vals)
}


# helper for KM method: allow log-scale mean calculation even if log-times shift < 0
calc_cenfit_negatives <- function(cf) {
  if (!inherits(cf, "cenfit")) stop("Input must be a cenfit object.")
  s <- cf@survfit
  if (!is.null(s$strata)) stop("Stratified cenfit objects not supported.")

  ev_times <- s$time[s$n.event > 0]
  if (length(ev_times) == 0) stop("No detected events to estimate from.")

  min_time <- min(ev_times, na.rm = TRUE)
  shift_amt <- if (min_time < 0) -min_time else 0

  cf2 <- if (shift_amt > 0) {
    s2 <- s
    s2$time <- s2$time + shift_amt
    new("cenfit", survfit = s2)
  } else {
    cf
  }

  m    <- mean(cf2)["mean"]
  sd_v <- stats::sd(cf2)

  if (shift_amt > 0) m <- m - shift_amt

  list(
    mean = unname(m),
    sd   = unname(sd_v)
  )
}


#' Kaplan–Meier estimate for left-censored data
#'
#' Uses `cenfit()` (bundled, NADA-derived) on the raw scale for mean/SD, and on
#' the log scale for geometric mean/GSD.
#'
#' @param censored_sim Numeric vector with observed values (LOD at censored indices).
#' @param censored_bool Logical vector (TRUE = censored).
#' @param ... Passed to [cenfit()].
#' @return Named list with mean, sd, gm, gsd.
#' @export
KM_estimate <- function(censored_sim, censored_bool, ...) {
  # KM cannot handle non-positive values for log-scale summary stats.
  if (any(censored_sim <= 0, na.rm = TRUE)) {
    return(list(mean = NA_real_, sd = NA_real_, gm = NA_real_, gsd = NA_real_))
  }

  mycenfit <- cenfit(obs = censored_sim, censored = censored_bool, ...)

  logcenfit <- cenfit(obs = log(censored_sim), censored = censored_bool, ...)
  accurate_log_vals <- calc_cenfit_negatives(logcenfit)

  list(
    mean = unname(mean(mycenfit)[1]),
    sd   = stats::sd(mycenfit),
    gm   = exp(accurate_log_vals[[1]]),
    gsd  = exp(accurate_log_vals[[2]])
  )
}


#' Simple substitution methods for left-censored data
#'
#' @param censored_sim Numeric vector with observed values (LOD at censored indices).
#' @param censored_bool Logical vector (TRUE = censored).
#' @param method One of `"root2"`, `"half"`, `"LOD"`, or `"zero"`.
#' @return Named list with mean, sd, gm, gsd.
#' @export
substitution_estimate <- function(censored_sim, censored_bool,
                                 method = c("root2", "half", "LOD", "zero")) {
  method <- match.arg(method)

  subs_vals <- censored_sim[censored_bool]

  censored_sim[censored_bool] <- switch(
    method,
    root2 = subs_vals / sqrt(2),
    half  = subs_vals / 2,
    LOD   = subs_vals,
    zero  = 0
  )

  estimate(censored_sim)
}


#' Fit a censored-data logspline density (Kooperberg–Stone 1992)
#'
#' This calls the modified [oldlogspline()] bundled with this package.
#'
#' @param censored_sim Numeric vector with observed values (LOD at censored indices).
#' @param censored_bool Logical vector (TRUE = censored).
#' @return An `oldlogspline` object, or `NULL` on failure.
#' @export
spline_fit <- function(censored_sim, censored_bool) {
  lbound <- 0
  ubound <- max(censored_sim, na.rm = TRUE)

  fit <- tryCatch(
    suppressMessages({
      oldlogspline(
        uncensored = censored_sim[!censored_bool],
        left       = censored_sim[censored_bool],
        lbound     = lbound,
        ubound     = ubound
      )
    }),
    error = function(e) {
      message("Logspline fit failed: ", e$message)
      NULL
    }
  )

  fit
}


#' Logspline moment estimation (Spline)
#'
#' Estimates the first two raw moments via numeric integration of the fitted
#' logspline density, then converts them to mean/SD. Also estimates geometric
#' mean/GSD by integrating on the log scale.
#'
#' @param censored_sim Numeric vector with observed values (LOD at censored indices).
#' @param censored_bool Logical vector (TRUE = censored).
#' @param fit An `oldlogspline` object, typically from [spline_fit()].
#' @return Named list with mean, sd, gm, gsd (may be `NA` on failure).
#' @export
spline_estimate <- function(censored_sim, censored_bool, fit) {
  # require vector length <14
  if (length(censored_sim) < 14) {
    return(list(mean = NA_real_, sd = NA_real_, gm = NA_real_, gsd = NA_real_))
  }

  if (any(censored_sim <= 0, na.rm = TRUE)) {
    warning("Non-positive values detected; returning NAs for all estimates.")
    return(list(mean = NA_real_, sd = NA_real_, gm = NA_real_, gsd = NA_real_))
  }

  if (is.null(fit)) {
    return(list(mean = NA_real_, sd = NA_real_, gm = NA_real_, gsd = NA_real_))
  }

  lower_bound <- 0
  upper_bound <- max(censored_sim, na.rm = TRUE)

  samp_first_moment  <- mean(censored_sim)
  samp_second_moment <- mean(censored_sim^2)

  pdf_x <- function(x) logspline::doldlogspline(x, fit)

  f_vec <- function(x) {
    if (is.matrix(x)) {
      xs    <- as.numeric(x[1, ])
      n_pts <- ncol(x)
    } else {
      xs    <- as.numeric(x)
      n_pts <- length(xs)
    }
    vals <- xs * pdf_x(xs)
    matrix(vals, nrow = 1, ncol = n_pts)
  }

  first_moment <- tryCatch({
    res1 <- cubature::hcubature(
      f               = f_vec,
      lowerLimit      = lower_bound,
      upperLimit      = upper_bound,
      vectorInterface = TRUE,
      maxEval         = 1e4,
      tol             = 1e-6,
      absError        = 1e-8
    )

    integ1 <- res1$integral

    if (!is.finite(integ1) || abs(integ1 - samp_first_moment) > 1e6) {
      warning("Integral invalid or differs by > 1e6; returning NA")
      NA_real_
    } else {
      integ1
    }
  }, error = function(e) {
    message("Error computing first moment: ", e$message)
    NA_real_
  })

  g_vec <- function(x) {
    if (is.matrix(x)) {
      xs    <- as.numeric(x[1, ])
      n_pts <- ncol(x)
    } else {
      xs    <- as.numeric(x)
      n_pts <- length(xs)
    }
    vals <- xs^2 * pdf_x(xs)
    matrix(vals, nrow = 1, ncol = n_pts)
  }

  second_moment <- tryCatch({
    res2 <- cubature::hcubature(
      f               = g_vec,
      lowerLimit      = lower_bound,
      upperLimit      = upper_bound,
      vectorInterface = TRUE,
      maxEval         = 1e4,
      tol             = 1e-6,
      absError        = 1e-8
    )

    integ <- res2$integral

    if (!is.finite(integ) || abs(integ - samp_second_moment) > 1e6) {
      warning("Integral invalid or differs by > 1e6; returning NA")
      NA_real_
    } else {
      integ
    }
  }, error = function(e) {
    message("Error computing second moment: ", e$message)
    NA_real_
  })

  var_orig <- second_moment - first_moment^2

  # ---------- log-scale moments via hcubature (vectorized) ----------
  log_upper <- log(upper_bound)
  log_lower <- -700

  samp_log_first_moment <- mean(log(censored_sim))
  samp_log_second_central_moment <- mean((log(censored_sim) - samp_log_first_moment)^2)

  l1_vec <- function(tmat) {
    if (is.matrix(tmat)) {
      ts   <- as.numeric(tmat[1, ])
      npts <- ncol(tmat)
    } else {
      ts   <- as.numeric(tmat)
      npts <- length(ts)
    }
    xs <- exp(ts)
    dens <- numeric(npts)
    keep <- xs > 0
    dens[keep] <- pdf_x(xs[keep])
    vals <- ts * dens * xs
    matrix(vals, nrow = 1, ncol = npts)
  }

  mu_log <- tryCatch({
    res_log1 <- cubature::hcubature(
      f               = l1_vec,
      lowerLimit      = log_lower,
      upperLimit      = log_upper,
      vectorInterface = TRUE,
      maxEval         = 1e4,
      tol             = 1e-6,
      absError        = 1e-8
    )

    integ_log1 <- res_log1$integral

    if (!is.finite(integ_log1) || abs(integ_log1 - samp_log_first_moment) > log(1e6)) {
      warning("Log first moment invalid or differs by > 1e6; returning NA")
      NA_real_
    } else {
      integ_log1
    }
  }, error = function(e) {
    message("Error computing mu_log (hcubature): ", e$message)
    NA_real_
  })

  l2c_vec <- function(tmat) {
    if (is.matrix(tmat)) {
      ts   <- as.numeric(tmat[1, ])
      npts <- ncol(tmat)
    } else {
      ts   <- as.numeric(tmat)
      npts <- length(ts)
    }
    xs <- exp(ts)
    dens <- numeric(npts)
    keep <- xs > 0
    dens[keep] <- pdf_x(xs[keep])
    vals <- (ts - mu_log)^2 * dens * xs
    matrix(vals, nrow = 1, ncol = npts)
  }

  var_log <- tryCatch({
    res_log2 <- cubature::hcubature(
      f               = l2c_vec,
      lowerLimit      = log_lower,
      upperLimit      = log_upper,
      vectorInterface = TRUE,
      maxEval         = 1e4,
      tol             = 1e-6,
      absError        = 1e-8
    )

    integ_log2 <- res_log2$integral

    if (!is.finite(integ_log2) || abs(integ_log2 - samp_log_second_central_moment) > log(1e6)) {
      warning("Log variance invalid or differs by > 1e6; returning NA")
      NA_real_
    } else {
      integ_log2
    }
  }, error = function(e) {
    message("Error computing var_log (hcubature): ", e$message)
    NA_real_
  })

  gm  <- if (is.na(mu_log)) NA_real_ else exp(mu_log)
  gsd <- if (is.na(var_log)) NA_real_ else exp(sqrt(var_log))

  list(
    mean = first_moment,
    sd   = sqrt(var_orig),
    gm   = gm,
    gsd  = gsd
  )
}


#' Logspline robust imputation estimator (rSpline)
#'
#' Computes Helsel–Cohn plotting positions and imputes censored values using
#' quantiles of the fitted censored-data logspline density.
#'
#' @param censored_sim Numeric vector with observed values (LOD at censored indices).
#' @param censored_bool Logical vector (TRUE = censored).
#' @param fit An `oldlogspline` object, typically from [spline_fit()].
#' @return Named list with mean, sd, gm, gsd (may be `NA` on failure).
#' @export
robust_spline <- function(censored_sim, censored_bool, fit) {
  if (length(censored_sim) != length(censored_bool)) {
    stop("'censored_sim' and 'censored_bool' must be the same length")
  }

  if (length(censored_sim) < 14) {
    return(list(mean = NA_real_, sd = NA_real_, gm = NA_real_, gsd = NA_real_))
  }

  if (is.null(fit)) {
    return(list(mean = NA_real_, sd = NA_real_, gm = NA_real_, gsd = NA_real_))
  }

  pp <- hc.ppoints(censored_sim, censored_bool)
  imputed_values <- logspline::qoldlogspline(pp[censored_bool], fit)

  all_vals <- censored_sim
  all_vals[censored_bool] <- imputed_values

  estimate(all_vals)
}


#' Beta substitution estimator (Hewett & Ganser)
#'
#' Implements the beta substitution approach described by Hewett and Ganser
#' (2007). Designed for lognormal-like concentration data with left-censoring.
#'
#' @param censored_sim Numeric vector with observed values (LOD at censored indices).
#' @param censored_bool Logical vector (TRUE = censored).
#' @param LOD Detection limit value(s). May be a scalar or vector.
#' @return Named list with mean, sd, gm, gsd.
#' @export
beta_sub <- function(censored_sim, censored_bool, LOD) {
  # calculate an average field LOD according to Hewett and Ganser
  if (length(LOD) > 1) {
    cens_vals <- censored_sim[censored_bool]
    mi <- vapply(LOD, function(dl) sum(cens_vals == dl), numeric(1))
    LOD <- exp(sum(mi * log(LOD)) / sum(mi))
  }

  n <- length(censored_sim)
  k <- sum(censored_bool)

  if (k == 0) {
    return(estimate(censored_sim))
  }
  if (k == n) {
    return(list(mean = NA_real_, sd = NA_real_, gm = NA_real_, gsd = NA_real_))
  }

  if (any(censored_sim <= 0)) stop("Values must be positive for log-transformation")

  ln_detects <- log(censored_sim[!censored_bool])

  y_hat <- (1 / (n - k)) * sum(ln_detects)
  z <- stats::qnorm((k / n))

  f_z <- stats::dnorm(z, 0, 1) / (1 - stats::pnorm(z, 0, 1))
  s_hat_y <- (y_hat - log(LOD)) / (f_z - z)

  f_sy_z <- (1 - stats::pnorm(z - (s_hat_y / n), 0, 1)) /
    (1 - stats::pnorm(z, 0, 1))

  beta_mean <- (n / k) * stats::pnorm(z - s_hat_y, 0, 1) *
    exp(-s_hat_y * z + (s_hat_y^2) / 2)

  sim_tmp <- censored_sim
  sim_tmp[censored_bool] <- beta_mean * LOD
  mean_est <- mean(sim_tmp)

  beta_GM <- exp(-(n - k) * n / k * log(f_sy_z) - s_hat_y * z - (n - k) / (2 * k * n) * s_hat_y^2)
  sim_tmp <- censored_sim
  sim_tmp[censored_bool] <- beta_GM * LOD
  gm_est <- exp(mean(log(sim_tmp)))

  ratio <- mean_est / gm_est

  if (is.na(ratio) || gm_est <= 0) {
    s_y <- NA_real_
  } else if (ratio <= 1) {
    s_y <- 0
  } else {
    s_y <- sqrt(2 * n / (n - 1) * log(ratio))
  }

  gsd_est <- exp(s_y)

  # Convert log-scale spread (s_y = sdlog) into arithmetic SD under a
  # lognormal assumption: SD = mean * CV, where CV = sqrt(exp(sdlog^2) - 1).
  sd_est <- if (is.na(s_y)) {
    NA_real_
  } else {
    mean_est * sqrt(exp(s_y^2) - 1)
  }

  list(
    mean = mean_est,
    sd   = sd_est,
    gm   = gm_est,
    gsd  = gsd_est
  )
}



