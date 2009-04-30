\name{regsimh}
\alias{regsimh}
\alias{print.regsimh}
\title{Simulate the distribution of heterogeneity and goodness-of-fit measures}
\description{
Estimates, using Monte Carlo simulation, the distribution of
heterogeneity and goodness-of-fit measures for regional frequency analysis.
These are the statistics \eqn{H} and \eqn{Z^{\rm DIST}}{Z^{DIST}}
defined respectively in sections 4.3.3 and 5.2.3 of Hosking and Wallis (1997).
}
\synopsis{
print.regsimh(x, ...)
}
\usage{
regsimh(qfunc, para, cor = 0, nrec, nrep = 500, nsim = 500)
}
\arguments{
  \item{qfunc}{List containing the quantile functions for each site.
        Can also be a single quantile function, which will be used for each site.}
  \item{para}{Parameters of the quantile functions at each site.
        If \code{qfunc} is a list, \code{para} must be a list of the same
        length whose components are numeric vectors, the parameters of
        the corresponding component of \code{qfunc}.
        If \code{qfunc} is a single quantile function, \code{para} can be
        a single vector, containing a single set of parameter values
        that will be used for each site;
        a matrix or data frame whose rows each contain the
        parameter values for one site;
        or a list of length \code{length(nrec)} whose components are
        numeric vectors, each containing the parameter values for one site.}
  \item{cor}{Specifies the correlation matrix of the frequency distribution
        of each site's data.  Can be a matrix (which will be rescaled to a
        correlation matrix if necessary) or a constant (which will be taken
        as the correlation between each pair of sites).}
  \item{nrec}{Numeric vector containing the record lengths at each site.}
  \item{nrep}{Number of simulated regions.}
  \item{nsim}{Number of simulations used, within each of the \code{nrep}
        simulated regions, when calculating heterogeneity and
        goodness-of-fit measures.}
}
\details{
A  realization is generated of data simulated from
the region specified by parameters \code{qfunc}, \code{para}, and \code{cor},
and with record lengths at each site specified by argument \code{nrec}.
The simulation procedure is as described in Hosking and Wallis (1997),
Table 6.1, through step 3.1.2.
Heterogeneity and goodness-of-fit measures are computed
for the realization, using the same method as in function \code{regtst}.
The entire procedure is repeated \code{nrep} times, and the values
of the heterogeneity and goodness-of-fit measures are saved.
Average values, across all \code{nrep} realizations,
of the heterogeneity and goodness-of-fit measures are computed.
}
\value{An object of class \code{"regsimh"}.
This is a list with the following components:
\item{nrep}{The number of simulated regions (argument \code{nrep}).}
\item{nsim}{The number of simulation used within each region
      (argument \code{nsim}).}
\item{results}{Matrix of dimension 8 \eqn{\times}{by} \code{nrep},
      containing the values, for each of the \code{nrep} simulated regions,
      of the heterogeneity and goodness-of-fit measures.}
\item{means}{Vector of length 8, containing the mean values,
      across the \code{nrep} simulated regions,
      of the three heterogeneity and five goodness-of-fit measures.}
}
\references{
  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
%\note{}
\seealso{\code{\link{regtst}} for details of the
heterogeneity and goodness-of-fit measures.}
\examples{
\dontrun{ % takes too much time
data(Cascades)            # A regional data set

rmom<-regavlmom(Cascades) # Regional average L-moments

# Set up an artificial region to be simulated:
# -- Same number of sites as Cascades
# -- Same record lengths as Cascades
# -- Mean 1 at every site (results do not depend on the site means)
# -- L-CV varies linearly across sites, with mean value equal
#    to the regional average L-CV for the Cascades data.
#    'LCVrange' specifies the  range of L-CV across the sites.
# -- L-skewness the same at each site, and equal to the regional
#    average L-skewness for the Cascades data
nsites <- nrow(Cascades)
means <- rep(1,nsites)
LCVrange <- 0.025
LCVs <- seq(rmom[2]-LCVrange/2, rmom[2]+LCVrange/2, len=nsites)
Lskews<-rep(rmom[3], nsites)

# Each site will have a generalized normal distribution:
# get the parameter values for each site
pp <- t(apply(cbind(means, means*LCVs ,Lskews), 1, pelgno))

# Set correlation between each pair of sites to 0.64, the
# average inter-site correlation for the Cascades data
avcor <- 0.64

# Run the simulation.  It will take some time (about 1 minute
# on a Lenovo T60, a moderately fast 2006-vintage laptop)
# Note that the results are consistent with the statement
# "the average H value of simulated regions is 1.08"
# in Hosking and Wallis (1997, p.98).
set.seed(123)
regsimh(qfunc=quagno, para=pp, cor=avcor, nrec=Cascades$n,
  nrep=100)
}}
\keyword{misc}
