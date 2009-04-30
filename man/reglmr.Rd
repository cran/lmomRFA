\name{reglmr}
\alias{reglmr}
\title{Regional weighted average of L-moments}
\description{
Regional weighted average of \eqn{L}-moments.
}
\usage{
reglmr(xmom, weight)
}
\arguments{
  \item{xmom}{
    Matrix or data frame each of whose rows contains the \eqn{L}-moments
    and \eqn{L}-moment ratios for one site, in the order \eqn{\ell_1}{l_1},
    \eqn{\ell_2}{l_2}, \eqn{t_3}{t_3}, \eqn{t_4}{t_4}, etc.
  }
  \item{weight}{
    Vector containing the weights to be used for each site.
    If omitted, equal weights will be used.
  }
}
%\details{}
\value{
Vector containing the regional average \eqn{L}-moments.
}
%\references{}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
\note{This function is deprecated and may be removed
from a future version of the package.
Function \code{\link{regavlmom}} is the recommended replacement.
}
%\seealso{}
\examples{
data(Maxwind)
(xmom<-t(sapply(Maxwind,samlmu)))
nrec<-sapply(Maxwind,length)
reglmr(xmom,nrec)   # weighted by record length
reglmr(xmom)        # unweighted
}
\keyword{misc}
