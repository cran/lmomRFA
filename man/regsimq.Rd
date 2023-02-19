\name{regsimq}
\alias{regsimq}
\alias{print.regsimq}
\title{Compute error bounds for a fitted regional frequency distribution}
\description{
Computes, using Monte Carlo simulation, relative error bounds
for estimated quantiles of a regional frequency distribution
fitted by \code{regfit}.
}
%\synopsis{
%print.regsimq(x, ...)
%}
\usage{
regsimq(qfunc, para, cor = 0, index = NULL, nrec, nrep = 10000,
        fit = "gev", f = c(0.01, 0.1, 0.5, 0.9, 0.99, 0.999),
        boundprob = c(0.05, 0.95), save = TRUE)
}
\arguments{
  \item{qfunc}{List containing the quantile functions for each site.
        Can also be a single quantile function, which will be used for each site.}
  \item{para}{Parameters of the quantile functions at each site.
        If \code{qfunc} is a list, \code{para} must be a list of the same length
        whose components are numeric vectors, the parameters of
        the corresponding component of \code{qfunc}.
        If \code{qfunc} is a single quantile function, \code{para} can be
        a single vector, containing a single set of parameter values
        that will be used for each site;
        a matrix or data frame whose rows each contain the
        parameter values for one site;
        or a list of length \code{length(nrec)} whose components are
        numeric vectors, each containing the parameter values for one site.}
  \item{cor}{Specifies the correlation matrix of the frequency distribution of
        each site's data.  Can be a matrix (which will be rescaled to a
        correlation matrix if necessary) or a constant (which will be taken as
        the correlation between each pair of sites).}
  \item{index}{Specifies the value of the site-specific scale factor
        (\dQuote{index flood}) at each site.  Can be:
        a vector, containing the values at each site;
        a constant, which will be taken to be the index flood value at each site;
        or (the default) \code{NULL}, in which case the index floods at each site
        will be taken to be the means of the quantile functions implied by
        \code{qfunc} and \code{para}, and will be computed
        by numerical integration of those quantile functions.}
  \item{nrec}{Numeric vector containing the record lengths at each site.}
  \item{nrep}{Number of simulated regions.}
  \item{fit}{Character string specifying the distribution to be fitted.
        See \dQuote{Details} below.}
  \item{f}{Vector of probabilities corresponding to the quantiles
        whose accuracy is to be estimated.}
  \item{boundprob}{Vector of probabilities for which error bounds will
        be computed.}
  \item{save}{Logical. Should the simulated values of the ratio of
        the estimated to the true regional growth curve be saved?
        These values are needed when \code{sitequantbounds} is called with
        its argument \code{index} present, e.g. to compute error bounds for
        quantiles at sites other than those whose data were used to fit the
        regional frequency distribution (e.g., ungauged sites).
        If this computation is not required, storage can be saved
        by setting \code{save} to \code{FALSE}.}
}
\details{
A  realization of data from a region is generated as follows.
The frequency distributions at sites (specified by
arguments \code{qfunc} and \code{para}) are expressed
as \eqn{Q_i(F)=\mu_i q_i(F)}{Q_i(F)=mu_i q_i(F)} where \eqn{\mu_i}{mu_i} is the site-specific scale factor
(\dQuote{index flood}) and \eqn{q_i(F)} is the at-site growth curve.
At each simulation run the at-site growth curves of each site
are randomly permuted, and are scaled by the (unpermuted)
index flood values for the sites.
Data are simulated from these frequency distributions,
with inter-site correlation specified by argument \code{cor}
and record lengths at each site specified by argument \code{nrec}.
The regional frequency distribution specified by argument \code{fit}
is then fitted to the simulated data, as in function \code{regfit}.
The procedure is as described in Hosking and Wallis (1997), Table 6.1,
except that the permutation of at-site growth curves is a later
modification, intended to give more realistic sets of simulated data.
For more details, including exact definitions of quantities computed
in the simulation and returned by functions \code{regsimq},
\code{regquantbounds}, and \code{regsitebounds}, see vignette \code{RegSim}.

From each realization the sample mean values and estimates of the
quantiles of the regional growth curve, for nonexceedance probabilities
specified by argument \code{f}, are saved.

From the simulated values, for each quantile specified by argument \code{f}
the relative root mean square error (relative RMSE) is computed as in
Hosking and Wallis (1997, eq. (6.15)).
Error bounds are also computed, as in Hosking and Wallis (1997, eq. (6.18))
but with bound probabilities specified by argument \code{boundprob}
rather than the fixed values 0.05 and 0.95 considered by Hosking and Wallis.
The error bounds are sample quantiles, across the \code{nrep} realizations,
of the ratio of the estimated regional growth curve
to the true at-site growth curve
or of the ratio of the estimated to the true quantiles at individual sites.

For distribution \code{fit} there should exist a function to estimate
the parameters of the distribution given a set of \eqn{L}-moments.
The function should have a name that is the character string
\code{"pel"} followed by the character string \code{fit}.
It should accept a single argument, a vector containing \eqn{L}-moments
\eqn{\ell_1}{l_1}, \eqn{\ell_2}{l_2}, \eqn{t_3}, \eqn{t_4}, etc.,
and return a vector of distribution parameters.

For distribution \code{fit} there should also exist a quantile function,
which should have a name that is the character string
\code{"qua"} followed by the character string \code{fit}.
It should accept two arguments: a vector of probabilities
and a vector containing the parameters of the distribution.

The search path used to find the \code{"pel"} and \code{"qua"} functions
is the same as for arguments supplied to \code{regsimq}, i.e.
the enclosing frames of the function, followed by the search path
specified by \code{search()}.

The estimation routines and quantile functions in package \pkg{lmom}
have the form described here.   For example, to use a
generalized extreme value distribution set \code{fit} to be
the string \code{"gev"}; then the fitting function \code{pelgev}
and the quantile function \code{quagev} will be used
(unless these functions have been masked by another object
on the search path).
}
\value{An object of class \code{"regsimq"}.
This is a list with the following components:
\item{f}{Vector of probabilities corresponding to the quantiles
      whose accuracy is to be estimated.  A copy of argument \code{f}.}
\item{boundprob}{Vector of probabilities corresponding to the error bounds.
      A copy of argument \code{boundprob}.}
\item{nrep}{Number of simulated regions.}
\item{relbounds.rgc}{Data frame containing the relative RMSE and
      error bounds for the regional growth curve. It has columns
      \code{"f"}, probabilities corresponding to each quantile,
      \code{"rel.RMSE"}, relative RMSE of quantiles  of regional growth curve,
      and, for each bound probability in \code{boundprob},
      a column giving the error bound (quantile of the empirical distribution
      of simulated values of the ratio of the estimated regional growth curve
      to the true at-site growth curve) for that bound probability.}
\item{relbounds.by.site}{List of \code{length(nrec)} data frames.
      Each data frame contains the relative RMSE and error bounds
      for quantiles at one site, and has the same structure
      as component \code{relbounds.rgc} of the return value.}
\item{true.asgc}{If \code{save} is \code{TRUE}, a matrix of dimension
      \code{length(f)} \eqn{\times}{by} \code{length(nrec)}, containing
      values of the at-site growth curves (quantile functions at each site,
      divided by the site-specific scale factors)
      for quantiles corresponding to probabilities in \code{f}.

      If \code{save} is \code{FALSE}, list element
      \code{true.asgc} is \code{NULL}.}
\item{sim.rgc}{If \code{save} is \code{TRUE}, a matrix of dimension
      \code{length(f)} \eqn{\times}{by} \code{nrep}, containing
      the simulated values of the estimated regional growth curve
      for quantiles corresponding to probabilities in \code{f}.

      If \code{save} is \code{FALSE}, list element
      \code{sim.rgc} is \code{NULL}.}
}
\references{
  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}
\note{
Error bounds for the regional growth curve apply to the
regional growth curve regarded as an estimator of the
growth curve for a randomly chosen site.
The growth curve for site \eqn{i} is defined to be
\eqn{q_i(\,.\,)=Q_i(\,.\,)/\mu_i}{q_i(.) = Q_i(.) / mu_i}
where \eqn{Q_i(\,.\,)}{Q_i(.)} is the site's quantile function
and \eqn{\mu_i}{mu_i} is the site-specific scale factor
(\dQuote{index flood}) for site \eqn{i}.
For each of the \code{nrep} simulated regions,
and each probability \eqn{F} in \code{f},
the regional growth curve \eqn{\hat{q}(F)}{qhat(F)} is estimated
and the ratios \eqn{\hat{q}(F)/q_i(F)}{qhat(F)/q_i(F)} are calculated.
The relative error bounds are empirical quantiles,
corresponding to the probabilities in \code{boundprob},
of the \code{nrep * length(nrec)} values of these ratios
obtained from the simulations.

This differs from the approach taken in Hosking and Wallis (1997),
Table 6.2 and Fig. 6.2, in which error bounds are computed regarding
the estimated regional growth curve as an estimator of the regional
average growth curve \eqn{q^{\rm R}(\,.\,)}{q^R(.)}, the harmonic mean
of the sites' growth curves (Hosking and Wallis, 1997, p. 102).

When the parent region is homogeneous, with identical frequency distributions
at each site (apart from a site-specific scale factor), the two
approaches give identical results.  For heterogeneous regions the
\dQuote{regard as estimator of randomly chosen site growth curve} approach
yields error bounds that are wider, but are more realistic given that the
ultimate aim of regional frequency analysis is estimation of quantiles
at individual sites.
}
\seealso{
\code{\link{regfit}}, for details of fitting a regional frequency distribution;
\code{\link{regquantbounds}} and \code{\link{sitequantbounds}}, for
converting the relative bounds returned by \code{regsimq} into absolute bounds
for quantiles of the regional growth curve or of the
frequency distributions at individual sites.}
\examples{
data(Cascades)              # A regional data set

rmom <- regavlmom(Cascades) # Regional average L-moments

# Fit generalized normal distribution to regional data
rfit <- regfit(Cascades, "gno")

# Set up an artificial region to be simulated:
# -- Same number of sites as Cascades
# -- Same record lengths as Cascades
# -- Same site means as Cascades
# -- L-CV varies linearly across sites, with mean value equal
#    to the regional average L-CV for the Cascades data.
#    'LCVrange' specifies the  range of L-CV across the sites.
# -- L-skewness the same at each site, and equal to the regional
#    average L-skewness for the Cascades data
nsites <- nrow(Cascades)
means <- Cascades$mean
LCVrange <- 0.025
LCVs <- seq(rmom[2]-LCVrange/2, rmom[2]+LCVrange/2, len=nsites)
Lskews<-rep(rmom[3], nsites)

# Each site will have a generalized normal distribution:
# get the parameter values for each site
pp <- t(apply(cbind(means, means*LCVs ,Lskews), 1, pelgno))

# Set correlation between each pair of sites to 0.64, the
# average inter-site correlation for the Cascades data
avcor <- 0.64

# Run the simulation.  To save time, use only 100 replications.
simq <- regsimq(qfunc=quagno, para=pp, cor=avcor, nrec=Cascades$n,
  nrep=100, fit="gno")

# Relative RMSE and error bounds for the regional growth curve
simq$relbounds.rgc

# Relative RMSE and error bounds for quantiles at site 3
simq$relbounds.by.site[[3]]
}
\keyword{misc}
