\name{Appalach}
\alias{Appalach}
\docType{data}
\title{Data for streamflow gaging stations in Appalachia}
\description{
Site characteristics and sample \eqn{L}-moments of annual maximum streamflow
for 104 gaging stations in Appalachia.
}
\usage{data(Appalach)}
\format{
  A data frame with 104 observations on the following 11 variables:
  \describe{
    \item{\code{siteid}}{Character vector: each site's Hydrologic Unit Code,
                         a unique identifier.}
    \item{\code{lat}}{Numeric vector: gage latitude, in degrees.}
    \item{\code{long}}{Numeric vector: gage longitude, in degrees west of the Greenwich Meridian.}
    \item{\code{area}}{Numeric vector: drainage basin area, in square miles.}
    \item{\code{elev}}{Numeric vector: gage elevation, in feet.}
    \item{\code{n}}{Numeric vector: record length.}
    \item{\code{mean}}{Numeric vector: sample mean.}
    \item{\code{t}}{Numeric vector: sample \eqn{L}-CV.}
    \item{\code{t_3}}{Numeric vector: sample \eqn{L}-skewness.}
    \item{\code{t_4}}{Numeric vector: sample \eqn{L}-kurtosis.}
    \item{\code{t_5}}{Numeric vector: sample \eqn{L}-moment ratio \eqn{t_5}.}
  }
}
\details{
The data in columns \code{lat}, \code{long}, \code{area} and \code{elev},
and the streamflow data used to compute the sample \eqn{L}-moments,
were obtained from \dQuote{Hydrodata} CD-ROMs (Hydrosphere, 1993),
which reproduce data from the U.S. Geological Survey's WATSTORE data files.
}
\source{
The file \code{appalach.dat} in the LMOMENTS Fortran package (Hosking, 1996).
}
\references{
Hydrosphere (1993).
\emph{Hydrodata CD-ROMs, vol. 4.0: USGS peak values}.
Hydrosphere Data Products, Boulder, Colo.

Hosking, J. R. M. (1996).
Fortran routines for use with the method of \eqn{L}-moments, Version 3.
Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.
}
\examples{
data(Appalach)
}
\keyword{datasets}
