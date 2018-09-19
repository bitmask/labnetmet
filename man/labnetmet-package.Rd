\name{labnetmet-package}
\alias{labnetmet-package}
\alias{labnetmet}
\docType{package}
\title{
\packageTitle{labnetmet}
}
\description{
\packageDescription{labnetmet}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{labnetmet}

\packageIndices{labnetmet}

Calculate distances with hamming_dist, subgraph_dist or trans_dist to calculate distances between graphs (supplied as adjacency matricies).  Hamming_dist will calculate the number of differences in direct edges between graphs.  Subgraph_dist will calculate the number of subgraphs that differ.  Trans_dist will calculate the number of paths that differ between the two graphs, such that the transitive closure of a graph is identical to the graph.  Plots can be generated directly with plot_dist.  

}
\author{
\packageAuthor{labnetmet}

Maintainer: \packageMaintainer{labnetmet}
}
\references{
}
\keyword{ package }
\seealso{
}
\examples{
a1 <- matrix( c(0, 1, 1, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(a1) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

a2 <- matrix( c(0, 1, 1, 1,  0, 0, 1, 0,  0, 0, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(a2) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

labnetmet::hamming_dist(a1, a2)
labnetmet::subgraph_dist(a1, a2)
labnetmet::trans_dist(a1, a2)

a3 <- matrix( c(0, 0, 1, 1,  0, 0, 1, 0,  0, 1, 0, 0,  0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)
dimnames(a3) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

a4 <- matrix( c(0, 1, 0, 1,  0, 0, 1, 0,  0, 1, 0, 0,  1, 1, 1, 1), nrow=4, ncol=4, byrow=TRUE)
dimnames(a4) <- list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

labnetmet::plot_dist(list(a1, a2, a3, a4), trans_dist, "~/filename.pdf")
}
