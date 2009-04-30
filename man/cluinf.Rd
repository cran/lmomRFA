\name{cluinf}
\alias{cluinf}
\title{Provide information about a hierarchical clustering}
\description{
Agglomerative hierarchical clustering procedures typically produce
a list of the clusters merged at each stage of the clustering.
\code{cluinf} uses this list to construct arrays that explicitly show
which cluster a given data point belongs to, and which data points belong
to a given cluster.
}
\usage{
cluinf(merge, nclust)
}
\arguments{
  \item{merge}{Matrix with 2 columns. The \eqn{i}th row contains
  the labels of the clusters merged at the \eqn{i}th merge.

  Can also be the object returned by a call to \code{\link{cluagg}}.}
  \item{nclust}{Number of clusters.}
}
%\details{}
\value{Information about the clustering that has \code{nclust} clusters.
It is a list with the following elements:
\item{assign}{Vector giving the assignment of items to clusters.}
\item{list}{List with \code{nclust} elements. Each element contains
            the labels of the items in one cluster.}
\item{num}{Vector of length \code{nclust}, containing the
           number of items in each cluster.}
}
\references{
  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{hosking@watson.ibm.com}}
%\note{}
\seealso{\code{\link{cluagg}}}
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
(cl<-cluagg(att))
# Details of the clustering with 7 clusters
cluinf(cl, 7)
}
\keyword{cluster}
