\name{clukm}
\alias{clukm}
\title{Cluster analysis via K-means algorithm}
\description{
Performs cluster analysis using the K-means algorithm.
}
\usage{
clukm(x, assign, maxit = 10, algorithm = "Hartigan-Wong")
}
\arguments{
  \item{x}{
    A numeric matrix (or a data frame with all numeric
    columns, which will be coerced to a matrix).  Contains the data:
    each row should contain the attributes for a single point.}
  \item{assign}{
    A vector whose distinct values indicate
    the initial clustering of the points.}
  \item{maxit}{Maximum number of iterations.}
  \item{algorithm}{
    Clustering algorithm.  Permitted values are the same as for
    \code{kmeans}.}
}
%\details{}
\value{An object of class \code{kmeans}.  For details see the help
for \code{\link{kmeans}}.}   % .... is this adequate?
\references{
  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
\note{
\code{clukm} is a wrapper for the \R function \code{kmeans}.
The only difference is that in \code{clukm} the user supplies an initial
assignment of sites to clusters (from which cluster centers are computed),
whereas in \code{kmeans} the user supplies the initial cluster centers
explicitly.
}
\seealso{\code{\link{kmeans}}}
\examples{
## Clustering of gaging stations in Appalachia, as in Hosking
## and Wallis (1997, sec. 9.2.3)
data(Appalach)
# Form attributes for clustering (Hosking and Wallis's Table 9.4)
att <- cbind(a1 = log(Appalach$area),
             a2 = sqrt(Appalach$elev),
             a3 = Appalach$lat,
             a4 = Appalach$long)
att <- apply(att, 2, function(x) x/sd(x))
att[,1] <- att[,1] * 3
# Clustering by Ward's method
(cl <- cluagg(att))
# Details of the clustering with 7 clusters
(inf <- cluinf(cl, 7))
# Refine the 7 clusters by K-means
clkm <- clukm(att, inf$assign)
# Compare the original and K-means clusters
table(Kmeans=clkm$cluster, Ward=inf$assign)
# Some details about the K-means clusters: range of area, number
# of sites, weighted average L-CV and L-skewness
bb <- by(Appalach, clkm$cluster, function(x)
  c( min.area = min(x$area),
     max.area = max(x$area),
     n = nrow(x),
     ave.t = round(weighted.mean(x$t, x$n), 3),
     ave.t_3 = round(weighted.mean(x$t_3, x$n), 3)))
# Order the clusters in increasing order of minimum area
ord <- order(sapply(bb, "[", 1))
# Make the result into a data frame.  Compare with Hosking
# and Wallis (1997), Table 9.5.
do.call(rbind, bb[ord])
}
\keyword{misc}
