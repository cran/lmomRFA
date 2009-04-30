\name{evplot.rfd}
\alias{evplot}
\alias{evplot.rfd}
\title{Extreme-value plot of a regional frequency distribution}
\description{
Plots a regional frequency distribution, optionally with error bounds
for either the regional growth curve
or the quantile function for an individual site.

The graph is an \dQuote{extreme-value plot}, i.e. the horizontal axis
is the quantile of an extreme-value type I (Gumbel) distribution,
and the quantile function of that distribution would plot as a straight line.
}
%\synopsis{}
\usage{
\method{evplot}{rfd}(y, ybounds, npoints=101, plim, xlim=c(-2,5), ylim,
  xlab=expression("Reduced variate,  " * -log(-log(italic(F)))),
  ylab="Quantile", rp.axis=TRUE, type="l", lty=c(1,2), col=c(1,1),
  ...)
}
\arguments{
  \item{y}{Object of class \code{\link{rfd}}, containing the specification
        of a regional frequency distribution.}
  \item{ybounds}{Optional. Object of class \code{rfdbounds} (typically created
        by \code{\link{regquantbounds}} or \code{\link{sitequantbounds}}),
        containing error bounds for quantile estimates for the
        regional frequency distribution specified by \code{y}.}
  \item{npoints}{Number of points to use in drawing the quantile function.
        The points are equally spaced along the x axis.}
  \item{plim}{X axis limits, specified as probabilities.}
  \item{xlim}{X axis limits, specified as values of the Gumbel reduced variate
        \eqn{-\log(-\log(F))}{-log(-log(F))}, where \eqn{F}{F} is the
        nonexceedance probability.
        Not used if \code{plim} is specified.}
  \item{ylim}{Y axis limits.}
  \item{xlab}{X axis label.}
  \item{ylab}{Y axis label.}
  \item{rp.axis}{Logical. Whether to draw the \dQuote{Return period} axis,
        a secondary horizontal axis.}
  \item{type}{Vector of plot types.  The first element is for the quantile
        function; subsequent elements are for the error bounds,
        and will be used cyclically until all lines are drawn.
        Interpreted in the same way as the \code{type} plotting parameter,
        i.e. \code{"l"} for lines, \code{"b"} for points connected by lines, etc.}
  \item{lty}{Vector of line types. The first element is for the quantile
        function; subsequent elements are for the error bounds,
        and will be used cyclically until all lines are drawn.}
  \item{col}{Vector of colors. The first element is for the quantile
        function; subsequent elements are for the error bounds,
        and will be used cyclically until all lines are drawn.}
  \item{\dots}{Additional parameters are passed to the plotting routine.}
}
\details{
If \code{ybounds} is missing, a graph is drawn of the quantile function
(regional growth curve) of the distribution specified by \code{y}.

If \code{ybounds} is present, it may contain error bounds for either
a regional growth curve or the quantile function at a single site.
This regional growth curve or site quantile function is plotted
using arguments \code{type[1]}, \code{lty[1]} and \code{col[1]}.
Then, in each case, error bounds are added to the plot.
The \code{ybounds} object typically contains,
for several probabilities specified by \code{ybounds$bounds},
error bounds corresponding to that probability
for several quantiles.
For the\code{j}th bound probability, the bounds for the various quantiles
will be joined by straight lines (so to obtain a smooth curve
there should be a lot of quantiles!), using graphics parameters
\code{type[j+1]}, \code{lty[j+1]} and \code{col[j+1]}.
}
%\value{}
%\references{}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
%\note{}
\seealso{\code{\link{regfit}}, which creates objects of class \code{"rfd"};
\code{\link{regquantbounds}} and \code{\link{sitequantbounds}},
which create objects of class \code{"rfdbounds"};
%\code{\link{evdistp}}, \code{\link{evdistq}}, and \code{\link{evpoints}},
\code{evdistp}, \code{evdistq}, and \code{evpoints},
all in package \pkg{lmom}, for adding further curves and points
to the plot.
}
\examples{
data(Cascades)                  # An object of class "regdata"
rfit <- regfit(Cascades, "gno") # Fit a generalized normal distribution

evplot(rfit)                    # Plot the regional growth curve

# Compute error bounds for quantile estimates.  We will
# (optimistically) generate bounds for a homogeneous region
# with the same frequency distribution as the one fitted to
# the Cascades data.
fval <- seq(.01, .99, by=.01)   # A lot of quantiles
simq <- regsimq(rfit$qfunc, nrec=Cascades$n, nrep=100, f=fval,
  fit=rfit$dist)

# Regional growth curve, and bounds
rbounds <- regquantbounds(simq, rfit)
evplot(rfit, rbounds)

# Quantile function for site 3, and bounds
sbounds <- sitequantbounds(simq, rfit, site=3)
evplot(rfit, sbounds)
}
\keyword{hplot}
