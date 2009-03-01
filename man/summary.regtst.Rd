\name{summary.regtst}
\alias{summary.regtst}
\alias{print.summary.regtst}
\title{Summary of test statistics for regional frequency analysis}
\description{
\code{summary} method for an object of class \code{"regtst"}.
}
\usage{
\method{summary}{regtst}(object,
  prob=c(0.01,0.02,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.98,0.99,0.999),
  conf=0.90, decimals=c(4,4,2,3), ...)

\method{print}{summary.regtst}(x, decimals, ...)
}
\arguments{
  \item{object}{
    An object of class \code{"regtst"}, usually the result of a call to
    \code{regtst}.}
  \item{x}{
    An object of class \code{"summary.regtst"}, usually the result of a
    call to \code{summary.regtst}.}
  \item{prob}{
    Nonexceedance probabilities for which quantile estimates
    should be printed.}
  \item{conf}{
    Confidence level for printing parameter and quantile estimates.
    These quantities will be printed only for distributions that
    give an adequate fit at the specified confidence level.}
  \item{decimals}{
    Vector of length 4.  The four elements specify the number of
    decimal places to be used when printing \eqn{L}-moment ratios,
    distribution parameters, test statistics,
    and quantile estimates, respectively.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\details{
  The printed output corresponds closely to that produced by
  function \code{REGTST} in the LMOMENTS Fortran package (Hosking, 1996).
}
\value{
\code{summary.regtst} and \code{print.summary.regtst} each return,
invisibly, an object of class \code{"summary.regtst"}, which is a list
with elements as for class \code{"regtst"}, plus the following elements:

\item{conf}{
  Confidence level --- the \code{conf} argument supplied to
  \code{summary.regtst}.}

\item{prob}{
  Vector of nonexceedance probabilities --- the \code{prob} argument
  supplied to \code{summary.regtst}.}

\item{quant}{
  Matrix with 6 rows and \code{length(prob)} columns, containing
  quantile estimates for the five candidate distributions and the
  Wakeby distribution.}

\item{decimals}{
  Vector of length 4.  Number of decimals to be used when printing
  an object of class \code{"summary.regtst"} if the \code{decimals}
  argument of \code{print.summary.regtst} is not specified.}

}
\references{
  Hosking, J. R. M. (1996).
  Fortran routines for use with the method of \eqn{L}-moments, Version 3.
  Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.
}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
%\note{}
\seealso{\code{\link{regtst}}}
\examples{
# An example from Hosking (1996).  Compare the output with
# the file 'cascades.out' in the LMOMENTS Fortran package at
# http://lib.stat.cmu.edu/lmoments/general (results will not
# be identical, because random-number generators are different).
data(Cascades)
summary(regtst(Cascades, nsim=500))
}
\keyword{misc}
