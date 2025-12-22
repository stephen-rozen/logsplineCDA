\docType{package}
\name{logsplineCDA-package}
\alias{logsplineCDA}
\alias{logsplineCDA-package}
\title{logsplineCDA: helpers for left-censored (LOD) concentration data}
\description{
Package \pkg{logsplineCDA} provides the core helper functions used in
Rozen and Canales (2025) \emph{A Logspline-Based Approach for Analyzing
Left-Censored Concentration Data}. The package provides estimation methods for 
environmental concentration data with left-censoring at one or more limits of 
detection (LOD).
}
\details{
The main exported functions implement:
\itemize{
  \item logspline-based moment estimation (\code{spline_estimate})
  \item logspline-based robust imputation (\code{robust_spline})
  \item common alternatives (ROS, MLE, KM, substitution, beta substitution)
}

A small subset of code derived from the CRAN package \pkg{NADA} is bundled
for internal use (ROS and survival-based helpers) without requiring
\pkg{NADA} as a dependency.
}
\author{
Stephen Rozen and Robert Canales (paper authors).
}
\references{
Rozen, S. and Canales, R. (2025). \emph{A Logspline-Based Approach for
Analyzing Left-Censored Concentration Data}.
}
\keyword{package}
