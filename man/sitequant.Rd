\name{sitequant}
\alias{siteqfunc}
\alias{sitequant}
\title{Quantiles and quantile functions for individual sites in a region}
\description{
Quantiles and quantile functions for individual sites in a region.
\code{sitequant} computes quantiles directly;
\code{siteqfunc} returns a function that will compute quantiles.
}
%\synopsis{}
\usage{
sitequant(f, rfd, sitenames, index, drop = TRUE)

siteqfunc(rfd, sitenames, index)
}
\arguments{
  \item{f}{Vector of probabilities.}
  \item{rfd}{Object of class \code{\link{rfd}}, containing the specification
    of a regional frequency distribution.}
  \item{sitenames}{Vector of site names.}
  \item{index}{Values of the site-specific scale factor (\dQuote{index flood})
    for the sites.}
  \item{drop}{Logical: if \code{TRUE} and there is only one site,
    or one probability value, the value returned from \code{sitequant}
    will be a vector rather than a matrix.}
}
\details{
If \code{index} and \code{sitenames} are both present,
they must have the same length, and will be taken to refer to sites
whose names are the elements of \code{sitename} and whose index-flood values
are the elements of \code{index}.

If \code{index} is present and \code{sitenames} is missing,
quantiles are computed for sites whose index-flood values
are the elements of \code{index}; if \code{index} has names, these
names will be used as the site names.

If \code{sitenames} is present and \code{index} is missing,
then quantiles will be computed for a subset of the
sites in the region specified by \code{rfd}.
\code{sitenames} will be used to select sites from the vector \code{rfd$index},
either by position or by name.

If \code{sitenames} and \code{index} are both missing,
then quantiles will be computed for all of the
sites in the region specified by \code{rfd}.
}
\value{
For \code{sitequant}, a matrix whose rows contain quantiles
for a single site, for the probabilities specified in \code{f}.
If \code{drop} is \code{TRUE} and the matrix has only one row or column,
it will be returned as a vector.

For \code{siteqfunc}, a function or a list of functions that each compute
quantiles for one site.
Each function takes a single argument, a vector of probabilities,
and returns a vector of quantiles.
}
%\references{}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
%\note{}
%\seealso{}
\examples{
data(Cascades)
rfit <- regfit(Cascades, "gno")   # Fit regional distribution

## Quantiles for:
# - sites in the Cascades data set, indexed by number
sitequant(c(0.9, 0.99, 0.999), rfit, sitenames=1:3)

# - sites in the Cascades data set, indexed by name
sitequant(c(0.9, 0.99, 0.999), rfit,
  sitenames=c("350304", "351433", "351862"))

# - other sites, with specified index floods
sitequant(c(0.9, 0.99, 0.999), rfit, index=c(80, 100))

# - other sites, with specified index floods and names
sitequant(c(0.9, 0.99, 0.999), rfit, index=c(80, 100),
  sitenames=c("Site 80", "Site 100"))

# - a single site, with drop=FALSE: result is a matrix
sitequant(c(0.9, 0.99, 0.999), rfit, sitenames=10, drop=FALSE)

# - a single site, with drop=TRUE (the default): result is a vector
sitequant(c(0.9, 0.99, 0.999), rfit, sitenames=10)

# Quantile function for site 10
qfunc10 <- siteqfunc(rfit, site=10)

# Compute quantiles via the quantile function
qfunc10(c(0.9, 0.99, 0.999))

# Plot the quantile function
evplot(qfunc=qfunc10)

}
\keyword{misc}
