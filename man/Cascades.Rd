\name{Cascades}
\docType{data}
\alias{Cascades}
\title{L-moments of annual precipitation totals}
\description{
\eqn{L}-moments of annual precipitation totals for the
\dQuote{North Cascades} region of Plantico et al. (1990).
}
\usage{data(Cascades)}
\format{
  An object of class \code{\link{regdata}}.
  It is a  data frame with 19 observations on the following 7 variables:
  \describe{
    \item{\code{name}}{Character vector: site identifier.}
    \item{\code{n}}{Numeric vector: record length.}
    \item{\code{mean}}{Numeric vector: sample mean.}
    \item{\code{t}}{Numeric vector: sample \eqn{L}-CV.}
    \item{\code{t_3}}{Numeric vector: sample \eqn{L}-skewness.}
    \item{\code{t_4}}{Numeric vector: sample \eqn{L}-kurtosis.}
    \item{\code{t_5}}{Numeric vector: sample \eqn{L}-moment ratio \eqn{t_5}.}
  }
}
\details{
The data are summary statistics of annual precipitation totals
at 19 sites in the northwest U.S.,
the \dQuote{North Cascades} region of Plantico et al. (1990).
The precipitation data were obtained from the Historical Climatology Network
(Karl et al., 1990).
}
\source{
The file \code{cascades.dat} in the LMOMENTS Fortran package (Hosking, 1996).
}
\references{
Hosking, J. R. M. (1996).
Fortran routines for use with the method of \eqn{L}-moments, Version 3.
Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.

Karl, T. R., Williams, C. N., Quinlan, F. T., and Boden, T. A. (1990).
United States historical climatology network (HCN) serial temperature
and precipitation data.
\emph{ORNL/CDIAC-30 NDP-019/R1}, Carbon Dioxide Information
Analysis Center, Oak Ridge National Laboratory, Oak Ridge, Tenn.

Plantico, M. S., Karl, T. R., Kukla, G., and Gavin, J. (1990).
Is recent climate change across the United States related to rising
levels of anthropogenic greenhouse gases?
\emph{Journal of Geophysical Research}, \bold{95}, 16617--16637.
}
\examples{
data(Cascades)
Cascades        # L-moment ratios, etc., for the Cascades data
lmrd(Cascades)  # L-moment ratio diagram for the Cascades data
}
\keyword{datasets}
