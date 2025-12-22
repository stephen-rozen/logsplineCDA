\name{lod_methods}
\alias{lod_methods}
\alias{ROS_estimate}
\alias{MLE_estimate}
\alias{KM_estimate}
\alias{substitution_estimate}
\alias{spline_fit}
\alias{spline_estimate}
\alias{robust_spline}
\alias{beta_sub}
\title{Estimation methods for left-censored (LOD) data}
\description{
Functions implementing the estimators compared in Rozen and Canales (2025).
All estimators assume left-censoring where censored observations are recorded
at their detection limit value(s).
}
\usage{
ROS_estimate(censored_sim, censored_bool)
MLE_estimate(censored_sim, censored_bool, method = c("mle", "bcmle", "rmle"))
KM_estimate(censored_sim, censored_bool, ...)
substitution_estimate(censored_sim, censored_bool,
  method = c("root2", "half", "LOD", "zero"))
spline_fit(censored_sim, censored_bool)
spline_estimate(censored_sim, censored_bool, fit)
robust_spline(censored_sim, censored_bool, fit)
beta_sub(censored_sim, censored_bool, LOD)
}
\arguments{
\item{censored_sim}{Numeric vector of observations; entries at censored indices are LOD values.}
\item{censored_bool}{Logical vector; TRUE indicates a censored observation.}
\item{method}{Estimator variant; see Details.}
\item{fit}{An \code{oldlogspline} object returned by \code{spline_fit()}.}
\item{LOD}{Detection limit value(s). May be a scalar or a vector.}
\item{...}{Additional arguments passed to \code{cenfit()} (KM estimator).}
}
\details{
\itemize{
  \item \code{ROS_estimate}: Regression on order statistics (ROS), using
    \code{ros()} (NADA-derived).
  \item \code{MLE_estimate}: Lognormal maximum likelihood estimation based on
    \code{EnvStats::elnormCensored(..., method = "mle")}, with options for a
    bias-corrected mean (bcMLE) and a robust imputation variant (rMLE).
  \item \code{KM_estimate}: Kaplan-Meier estimator using \code{cenfit()}
    (NADA-derived).
  \item \code{substitution_estimate}: Simple substitutions (LOD, LOD/2,
    LOD/sqrt(2), or zero).
  \item \code{spline_fit}: Generates an \code{oldlogspline} object for use with 
  \code{spline_estimate} or \code{robust_spline}
  \item \code{spline_estimate}: Logspline-based moment
    estimation (\emph{Spline}) from a censored-data logspline density.
  \item \code{robust_spline}: Logspline-based robust imputation (\emph{rSpline})
    using logspline quantiles.
  \item \code{beta_sub}: Beta substitution estimator (Hewett \& Ganser).
}
}
\value{
Each estimator returns a named list with elements \code{mean}, \code{sd},
\code{gm}, and \code{gsd} (with \code{gsd} sometimes returned as \code{NA}
for methods that do not estimate it).
}

\examples{
# Example 1: Single-LOD data and the logspline moment estimator (Spline)
set.seed(1)
x <- stats::rlnorm(60, meanlog = log(2), sdlog = log(2))
LOD <- as.numeric(stats::quantile(x, probs = 0.25))

# Represent left-censoring by replacing non-detects with the LOD
censored_bool <- x < LOD
censored_sim <- x
censored_sim[censored_bool] <- LOD

fit <- spline_fit(censored_sim, censored_bool)
if (!is.null(fit)) {
  # Numerical integration over the fitted logspline density
  spline_estimate(censored_sim, censored_bool, fit)
}

# Example 2: Compare bcMLE and beta-substitution on the same censored dataset
MLE_estimate(censored_sim, censored_bool, method = "bcmle")
beta_sub(censored_sim, censored_bool, LOD = LOD)
}
\keyword{censored}
\keyword{distribution}
