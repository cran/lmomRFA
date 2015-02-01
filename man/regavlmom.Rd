\name{regavlmom}
\alias{regavlmom}
\title{Regional weighted average of L-moments}
\description{
Computes a regional weighted average of \eqn{L}-moments.
}
\usage{
  regavlmom(regdata, weight)
}
\arguments{
  \item{regdata}{Object of class \code{\link{regdata}}, containing
        summary statistics of the data for the sites in a region.}
  \item{weight}{Vector containing the weights to be used for each site.
        If omitted, weights will be the sample size at each site,
        taken from the second column of \code{regdata}.
        If a single value, equal weights will be used.}
}
%\details{}
\value{
Vector containing the regional average \eqn{L}-moments.
}
%\references{}
\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}
%\note{}
%\seealso{}
\examples{
reglm <- regsamlmu(Maxwind)
regavlmom(reglm)            # Weight proportional to record length
regavlmom(reglm, weight=1)  # Equal weights
}
\keyword{misc}
