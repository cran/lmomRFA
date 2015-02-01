\name{regquant}
\alias{regqfunc}
\alias{regquant}
\title{Quantiles and quantile function of a regional frequency distribution}
\description{
\code{regquant} computes quantiles of a regional frequency distribution,
i.e., values of the regional growth curve.

\code{regqfunc} returns a function that will compute the quantiles.
}
\usage{
regquant(f, rfd)

regqfunc(rfd)
}
\arguments{
  \item{f}{Vector of probabilities.}
  \item{rfd}{Object of class \code{\link{rfd}}, containing the specification
             of a regional frequency distribution.}
}
%\details{}
\value{\code{regquant} returns a vector of quantiles.

\code{regqfunc} returns the \code{qfunc} element of \code{rfd}.
This is a function that takes one argument,
which should be a vector of probabilities, and returns a vector of quantiles.
}
%\references{}
\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}
%\note{}
\seealso{\code{\link{regfit}}.}
\examples{
rfit <- regfit(Cascades,"gno")  # Fit regional distribution

# Compute some quantiles
regquant(seq(0.1, 0.9, by=0.1), regfit(Cascades,"gno"))

# Get the quantile function (regional growth curve)
rgc <- regqfunc(rfit)

# Compute quantiles by evaluating the regional growth curve
rgc(seq(0.1, 0.9, by=0.1))

# Plot the regional growth curve
curve(rgc, 0.01, 0.99)

}
\keyword{misc}
