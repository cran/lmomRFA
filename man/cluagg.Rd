\name{cluagg}
\alias{cluagg}
\title{Hierarchical clustering}
\description{
Performs cluster analysis by one of several agglomerative hierarchical methods.
}
\usage{
cluagg(x, method="ward")
}
\arguments{
  \item{x}{
    A numeric matrix (or a data frame with all numeric columns,
    which will be coerced to a matrix).  Contains the data:
    each row should contain the attributes for a single point.}
  \item{method}{Clustering method. Any method valid for
    \code{\link{hclust}} may be used.}
}
\details{
In agglomerative hierarchical clustering, there are initially \eqn{n} clusters,
each containing one data point, labeled \eqn{1} through \eqn{n}
in the same order as the data points.
At each stage of clustering, two clusters are merged.
Their labels are saved in the \code{merge} array.
The smaller of the two labels is used as the label of the merged cluster.
After the \eqn{i}th stage of clustering there are \eqn{n-i} clusters.
To find which data points belong to which clusters, use function
\code{\link{cluinf}}.
}
\value{A list with elements as follows.
\item{merge}{Matrix of dimension \code{(nrow(x)-1,2)}.
The \eqn{i}th row contains the labels of the clusters merged at the
\eqn{i}th merge.}
\item{wgss}{Vector of length \code{nrow(x)-1}.
The \eqn{i}th element is the total within-cluster dispersion
after the \eqn{i}th merge.}
}
\references{
  Hosking, J. R. M., and Wallis, J. R. (1997).
  \emph{Regional frequency analysis: an approach based on \eqn{L}-moments}.
  Cambridge University Press.
}
\author{J. R. M. Hosking \email{jrmhosking@gmail.com}}
\note{Clustering is performed internally by function \code{\link{hclust}}
in the \R \pkg{stats} package.}
\seealso{\code{\link{cluinf}} to get details of the clusters
at a particular stage of the merging.}
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
cluinf(cl,7)
}
\keyword{cluster}
