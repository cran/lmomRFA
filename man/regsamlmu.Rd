\name{regsamlmu}
\alias{regsamlmu}
\title{Sample L-moments of multiple data sets}
\description{
Computes the \dQuote{unbiased} sample \eqn{L}-moments and \eqn{L}-moment
ratios of multiple sets of data stored in a list or matrix.
Following the paradigm of regional frequency analysis,
we regard the data sets as coming from different measurement sites.
}
\usage{
regsamlmu(x, nmom = 5, sort.data = TRUE, lcv = TRUE)
}
\arguments{
  \item{x}{A list of numeric vectors, or a numeric matrix.}
  \item{nmom}{Number of \eqn{L}-moments to be computed.}
  \item{sort.data}{Whether each data set should be sorted.}
  \item{lcv}{Logical.
    If \code{TRUE}, the second \eqn{L}-moment will be expressed
    as a fraction of the mean, i.e. the computed value will be
    the sample \eqn{L}-CV \eqn{t=\ell_2/\ell_1}{t=l_2/l_1}.
    If \code{FALSE}, the second \eqn{L}-moment will simply be
    the sample \eqn{L}-scale value \eqn{\ell_2}{l_2}.}
}
\details{
  Sample \eqn{L}-moments are computed for each data set.
  The calculations use \code{samlmu} internally.
  If \code{x} is a list, each list element should contain data for one site
  and the names of the list elements should be the site names.
  If \code{x} is a matrix, each column should contain data for one site
  and the column names should be the site names.
}
\value{
  An object of class \code{\link{regdata}}.
  It is a data frame with columns \code{"name"} and \code{"n"},
  containing respectively the site names and the
  number of non-missing data values at each site,
  and further columns containing the \eqn{L}-moments and \eqn{L}-moment ratios,
  in the order \eqn{\ell_1}{l_1}, \eqn{t} (or \eqn{\ell_2}{l_2}),
  \eqn{t_3}{t_3}, \eqn{t_4}{t_4}, etc.
}

\note{The default parameter values are chosen to be convenient
for the regional frequency analysis methods described by
Hosking and Wallis (1997).
Note that the number of \eqn{L}-moments and the choice
of whether to return \eqn{L}-CV or \eqn{L}-scale
are different from the defaults for \code{samlmu}.
}

\references{
  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}

\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}

\examples{
data(Maxwind)       # a list
regsamlmu(Maxwind)

data(airquality)    # a data frame
regsamlmu(airquality[1:4])
}
\keyword{misc}
