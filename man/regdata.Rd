\name{regdata}
\alias{regdata}
\alias{as.regdata}
\title{The regdata class}
\description{
An object of class \code{"regdata"} stores
summary statistics of the data for the sites in a region.
It is a data frame with each row containing data for one site.
The columns should contain
the site name, record length and \eqn{L}-moments
and \eqn{L}-moment ratios, in the order
\eqn{\ell_1}{l_1} (mean),
\eqn{t} (\eqn{L}-CV),
\eqn{t_3}{t_3} (\eqn{L}-skewness),
\eqn{t_4}{t_4} (\eqn{L}-kurtosis),
\eqn{t_5}{t_5},
\eqn{t_6}{t_6}, etc.

There should be at least four columns, but most functions that use
objects of class \code{"regdata"} typically require more columns.
Six or seven columns (4 or 5 \eqn{L}-moments) is usually
adequate for regional frequency analysis.

Note that the fourth column should contain values of
the \eqn{L}-CV \eqn{t}, not the \eqn{L}-scale \eqn{\ell_2}{l_2}!

Objects of class \code{"regdata"} are created by \code{as.regdata},
and by \code{regsamlmu} (with default settings of its arguments).
They are used by several functions in package \pkg{lmomRFA}, including
\code{\link{regavlmom}} (which computes regional average \eqn{L}-moments),
\code{\link{regfit}} (which fits a regional frequency distribution), and
\code{\link{regtst}} (which computes discordancy, heterogeneity and
goodness-of-fit measures).
}
\usage{
as.regdata(x, warn.names=TRUE)
}
\arguments{
  \item{x}{\R object.}
  \item{warn.names}{Logical: if \code{TRUE}, warnings are issued if the
        column names of \code{x} appear to be inconsistent with
        what is expected for an object of class \code{"regdata"}.}
}
\details{
\code{as.regdata} converts an \R object to class \code{"regdata"}.
Only data frames and numeric matrices can be converted.
}
\value{An object of class \code{"regdata"}.}
%\references{}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
%\note{}
%\seealso{}
\examples{
data(Cascades)  # 'Cascades' is of class "regdata"
Cascades

# Create a data frame with site statistics
dd<-data.frame(
  name      =c("site 1", "site 2", "site 3"),
  n         =c(  20,   30,   40),
  mean      =c( 100,  110,  120),
  LCV       =c(0.20, 0.25, 0.30),
  L_skewness=c(0.15, 0.20, 0.25),
  L_kurtosis=c(0.10, 0.15, 0.20),
  t_5       =c(0.10, 0.12, 0.14))
# Convert to class "regdata"
rdd<-as.regdata(dd)
rdd
class(rdd)
}
\keyword{misc}
