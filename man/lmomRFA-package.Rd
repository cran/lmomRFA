\name{lmomRFA-package}
\docType{package}
\alias{lmomRFA}
\alias{lmomRFA-package}
\title{The lmomRFA package}
\description{
  \R functions for regional frequency analysis using \eqn{L}-moments.
}
\details{
This package implements methods described in the book
\dQuote{Regional frequency analysis: an approach based on \eqn{L}-moments}
by J. R. M. Hosking and J. R. Wallis.
It is a supplement to the \pkg{lmom} package,
which implements \eqn{L}-moment methods for more
general statistical applications.

The following functions are contained in this package.

\code{\link{cluagg}}
  performs agglomerative hierarchical clustering.

\code{\link{cluinf}}
  provides information about cluster membership in a
  hierarchical clustering.

\code{\link{clukm}}
  performs cluster analysis via the K-means algorithm.

\code{\link{as.regdata}}
  creates an object of class \code{"regdata"}, which contains a
  \dQuote{regional data set} consisting of summary statistics
  for different data samples, one of the main building blocks
  of regional frequency analysis.

\code{\link{regsamlmu}}
  computes the sample \eqn{L}-moments of multiple data sets.

\code{\link{regavlmom}} and \code{\link{reglmr}}
  both compute, with slightly different interfaces, a regional
  weighted average of sample \eqn{L}-moments from multiple sites.
  Function \code{regavlmom} is recommended for general use;
  \code{reglmr} is deprecated.

\code{\link{regtst}}
  computes discordancy, heterogeneity and goodness-of-fit measures
  for regional frequency analysis.  These statistics are as described
  in Hosking and Wallis (1997, chaps. 3-5).

\code{\link{regfit}}
  fits a frequency distribution to a regional data set, giving a
  \dQuote{regional frequency distribution}.

\code{\link{regqfunc}} and \code{\link{siteqfunc}}
  return the regional growth curve and the quantile functions for individual
  sites, respectively, from a regional frequency distribution fitted by
  \code{regfit}.

\code{\link{regquant}} and \code{\link{sitequant}}
  directly compute quantiles of the regional growth curve and of
  distributions for individual sites, respectively,
  from a regional frequency distribution fitted by \code{regfit}.

\code{\link{regsimh}}
  runs Monte Carlo simulations to estimate the distribution of
  heterogeneity and goodness-of-fit measures for an artificial region.

\code{\link{regsimq}}
  runs Monte Carlo simulations to estimate the variability of
  quantile estimates from a regional frequency distribution.

\code{\link{regquantbounds}} and \code{\link{sitequantbounds}}
  compute error bounds for the regional growth curve and for quantiles
  at individual sites, respectively, from a regional frequency distribution
  fitted by \code{regfit}.

Functions \code{cluagg}, \code{cluinf}, \code{clukm}, \code{reglmr},
and \code{regtst} are analogous to Fortran routines from the
LMOMENTS package, version 3.04, available from StatLib at
\url{https://lib.stat.cmu.edu/general/lmoments}.
In addition, functions \code{regsimh} and \code{regsimq}
provide similar functionality to \code{PROGRAM XSIM}
in the LMOMENTS Fortran package.
}
\references{
  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}
\keyword{package}
