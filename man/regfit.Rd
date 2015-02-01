\name{regfit}
\alias{regfit}
\alias{rfd}
\alias{print.rfd}
\title{Fit a regional frequency distribution}
\description{
Fits a frequency distribution to a vector of regional average
\eqn{L}-moments.  Returns an object of class \code{"rfd"},
which contains the specification of the regional frequency distribution:
the quantile function, parameters of the regional growth curve,
and the index-flood values (site-specific scale factors) for each site.
}
%\synopsis{
%print.rfd(x, ...)
%}
\usage{
regfit(regdata, dist)
}
\arguments{
  \item{regdata}{Object of class \code{\link{regdata}}, containing
        summary statistics of the data for the sites in a region.}
  \item{dist}{Character string specifying the distribution to be fitted.
        See \dQuote{Details} below.}
}
\details{
The function computes regional average \eqn{L}-moments
(by calling \code{regavlmom}) and fits a probability distribution
to the regional average \eqn{L}-moments.

For distribution \code{dist} there should exist a function to estimate
the parameters of the distribution given a set of \eqn{L}-moments.
The function should have a name that is the character string \code{"pel"}
followed by the character string \code{dist}.
It should accept a single argument, a vector containing \eqn{L}-moments
\eqn{\ell_1}{l_1}, \eqn{\ell_2}{l_2}, \eqn{t_3}, \eqn{t_4}, etc.,
and return a vector of distribution parameters.

For distribution \code{dist} there should also exist a quantile function,
which should have a name that is the character string
\code{"qua"} followed by the character string \code{dist}.
It should accept two arguments: a vector of probabilities
and a vector containing the parameters of the distribution.

The search path used to find the \code{"pel"} and \code{"qua"} functions
is the same as for arguments supplied to \code{regfit}, i.e.
the enclosing frames of the function, followed by the search path
specified by \code{search()}.

The estimation routines and quantile functions in package \pkg{lmom}
have the form described here.   For example, to use a
generalized extreme value distribution set \code{dist} to be
the string \code{"gev"}; then the fitting function \code{pelgev}
and the quantile function \code{quagev} will be used
(unless these functions have been masked by another object
on the search path).
}
\value{An object of class \code{"rfd"},
containing the specification of the regional frequency distribution:
It is a list with the following elements:
\item{dist}{The character string \code{dist}.}
\item{para}{Vector containing the parameters of the fitted regional distribution.}
\item{qfunc}{The quantile function of distribution \code{dist}.
It is a function that takes a single argument, a vector of probabilities,
and returns a vector of quantiles.}
\item{rmom}{The regional average \eqn{L}-moments.}
\item{index}{Index flood values at each site.  This is a named vector
whose values are the index flood values at each site, from \code{regdata[[3]]},
and whose names are the site names, from \code{regdata[[1]]}.}
}
%\references{}
\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}
%\note{}
%\seealso{}
\examples{
data(Cascades)                  # An object of class "regdata"
rfit <- regfit(Cascades, "gno") # Fit a generalized normal distribution
rfit                            # Print details of the fitted distribution
                                #   (components 'dist' and 'para')
rfit$index                      # Index flood values

evplot(rfit)                    # Plot the regional growth curve
evplot(qfunc=rfit$qfunc)        # The same, but with more typing
evplot(qfunc=regqfunc(rfit))    # The same, with still more typing
}
\keyword{misc}
