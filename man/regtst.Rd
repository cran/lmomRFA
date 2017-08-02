\name{regtst}
\alias{regtst}
\alias{regtst.s}
\alias{print.regtst}
\title{Test statistics for regional frequency analysis}
\description{
Computes discordancy, heterogeneity and goodness-of-fit measures
for regional frequency analysis.
These are the statistics \eqn{D_i}{D_i}, \eqn{H}, and \eqn{Z^{\rm DIST}}{Z^{DIST}}
defined respectively in sections 3.2.3, 4.3.3, and 5.2.3 of
Hosking and Wallis (1997).
}
\usage{
regtst(regdata, nsim=1000)

regtst.s(regdata, nsim=1000)
}
\arguments{
  \item{regdata}{
    Object of class \code{\link{regdata}} containing the input data.
    It should be a data frame, each of whose rows contains data for one site.
    The first seven columns should contain respectively
    the site name, record length and \eqn{L}-moments
    and \eqn{L}-moment ratios, in the order
    \eqn{\ell_1}{l_1} (mean),
    \eqn{t} (\eqn{L}-CV),
    \eqn{t_3}{t_3} (\eqn{L}-skewness),
    \eqn{t_4}{t_4} (\eqn{L}-kurtosis),
    and \eqn{t_5}{t_5}.

    Note that the fourth column should contain values of
    the \eqn{L}-CV \eqn{t}, not the \eqn{L}-scale \eqn{\ell_2}{l_2}!

    Function \code{regsamlmu}, with default settings
    of its arguments, returns an object of class \code{"regdata"}.}
  \item{nsim}{
    Number of simulations to use in the calculation of the
    heterogeneity and goodness-of-fit measures.

    If less than 2, only the discordancy measure will be calculated.
    }
}
\details{
The discordancy measure \eqn{D_i} indicates, for site \eqn{i},
the discordancy between the site's \eqn{L}-moment ratios
and the (unweighted) regional average \eqn{L}-moment ratios.
Large values might be used as a flag to indicate potential errors
in the data at the site.  \dQuote{Large} might be 3 for regions with 15
or more sites, but less (exact values in list element \code{Dcrit})
for smaller regions.

Three heterogeneity measures are calculated, each based on
a different measure of between-site dispersion of \eqn{L}-moment ratios:
[1] weighted standard deviation of \eqn{L}-CVs;
[2] average of \eqn{L}-CV/\eqn{L}-skew distances;
[3] average of \eqn{L}-skew/\eqn{L}-kurtosis distances.
These dispersion measures are the quantities \eqn{V}, \eqn{V_2},
and \eqn{V_3} defined respectively in equations (4.4), (4.6), and (4.7)
of Hosking and Wallis (1997).
The heterogeneity measures are calculated from them as in
equation (4.5) of Hosking and Wallis (1997).
In practice \code{H[1]} is probably sufficient.  A value greater than
(say) 1.0 suggests that further subdivision of the region should
be considered as it might improve the accuracy of quantile estimates.

Goodness of fit is evaluated for five candidate distributions:
generalized logistic,
generalized extreme value,
generalized normal (lognormal),
Pearson type III (3-parameter gamma), and
generalized Pareto.
In the output the distributions are referred to by 3-letter abbreviations,
respectively \code{glo}, \code{gev}, \code{gno}, \code{pe3}, and \code{gpa}.
If the region is homogeneous and data at different sites are
statistically independent, then if one of the distributions is
the true distribution for the region its goodness-of-fit measure
should have approximately a standard normal distribution.
Provided that the region is acceptably close to homogeneous,
the fit may be judged acceptable at the 10 per cent significance level
if the \eqn{Z} value is less than 1.645 (i.e., \code{qnorm(0.95)}) in absolute value.

Calculation of heterogeneity and goodness-of-fit measures
involves the sampling variability of \eqn{L}-moment ratios
in a homogeneous region whose record lengths and
average \eqn{L}-moment ratios match those of the data.
The sampling variability is estimated by Monte Carlo simulation
using \code{nsim} replications of the region.
Results will vary between invocations of \code{regtst}
with different seeds for the random-number generator.

In the homogeneous region used in the simulations, the sites have a
kappa distribution, fitted to the regional average \eqn{L}-moment ratios
of the data in \code{regdata}.  The kappa fit may fail if the regional average
\eqn{L}-kurtosis is high relative to the regional average \eqn{L}-skewness.
In this case a kappa distribution is fitted with shape parameter
\eqn{h} constrained to be \eqn{-1} (i.e., a generalized logistic distribution);
this gives the largest possible \eqn{L}-kurtosis value for a kappa distribution
with given \eqn{L}-skewness.

\code{regtst} and \code{regtst.s} are functionally identical.
\code{regtst} calls a Fortran routine internally and is faster,
typically by a factor of 3 or 4.
\code{regtst.s} is written almost entirely in the S language;
it is provided so that users can see how the calculations are done,
and can conveniently alter the code for their own purposes if necessary.
}
\value{
An object of class \code{"regtst"}, which is a list with elements as follows.

\item{data}{The input data, i.e. data frame \code{regdata}
  after coercion to class \code{"regdata"} if necessary.}
\item{nsim}{Number of simulations, i.e. the argument \code{nsim}.}
\item{D}{Vector containing the discordancy measures for each site.}
\item{Dcrit}{
  Vector of length 2 containing critical values of the discordancy
  measure corresponding to significance levels of 10 and 5 per cent
  --- except that the values never exceed 3 and 4 respectively.
  See Hosking and Wallis (1997), section 3.2.4.}
\item{rmom}{
  Vector of length 5 containing the regional weighted average
  \eqn{L}-moment ratios (weights proportional to record lengths).}
\item{rpara}{
  Vector of length 4 containing the parameters of a kappa distribution
  fitted to the regional weighted average \eqn{L}-moment ratios.}
\item{vobs}{
  Vector of length 3 containing the observed values of the three
  measures of between-site dispersion of \eqn{L}-moment ratios.}
\item{vbar}{
  Vector of length 3 containing the mean of the simulated values
  of the three dispersion measures.}
\item{vsd}{
  Vector of length 3 containing the standard deviation of the
  simulated values of the three dispersion measures.}
\item{H}{
  Vector of length 3 containing the three measures of regional
  heterogeneity.}
\item{para}{
  List of length 6 containing the parameters of the five candidate
  distributions and the Wakeby distribution (3-letter abbreviation
  \code{"wak"}) fitted to the regional weighted average
  \eqn{L}-moment ratios.}
\item{t4fit}{
  Vector of length 5 containing the \eqn{L}-kurtosis of the
  five candidate distributions fitted to the regional
  weighted average \eqn{L}-moment ratios.}
\item{Z}{
  Vector of length 5 containing the goodness-of-fit measures
  for each of the five candidate distributions.}
}
\references{
  Hosking, J. R. M. (1996).
  Fortran routines for use with the method of \eqn{L}-moments, Version 3.
  Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.

  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}
\note{
  Data frame \code{regdata} may have only six columns,
  i.e. the fifth \eqn{L}-moment ratio \eqn{t_5} may be omitted.
  In this case the return value will contain missing values for
  \code{rmom[5]} and the elements of \code{para$wak}.
}
\seealso{
\code{summary.regtst} for summaries.
}
\examples{
# An example from Hosking (1996).  Compare the output with
# the file 'cascades.out' in the LMOMENTS Fortran package at
# http://lib.stat.cmu.edu/lmoments/general (results will not
# be identical, because random-number generators are different).
summary(regtst(Cascades, nsim=500))

# Output from 'regsamlmu' can be fed straight into 'regtst'
regtst(regsamlmu(Maxwind))
}
\keyword{misc}
