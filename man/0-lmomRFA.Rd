\name{lmomRFA}
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

\code{\link{reglmr}}
  computes a regional weighted average of \eqn{L}-moments.

\code{\link{regsamlmu}}
  computes the sample \eqn{L}-moments of multiple data sets.

\code{\link{regtst}}
  computes discordancy, heterogeneity and goodness-of-fit measures
  for regional frequency analysis.  These statistics are as described
  in Hosking and Wallis (1997, chaps. 3-5).

Functions \code{cluagg}, \code{cluinf}, \code{clukm}, \code{reglmr},
and \code{regtst} are analogous to Fortran routines from the
LMOMENTS package, version 3.04, available from StatLib at
\url{http://lib.stat.cmu.edu/general/lmoments}.
}
\references{
  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
\keyword{univar}
