\name{plot_joint_graph}
\alias{plot_joint_graph}
\title{Plot joint ClustIRR graph}
\usage{plot_joint_graph(clust_irrs, cores = 1, as_visnet = FALSE)}
\arguments{
\item{clust_irrs}{list of two or mroe S4 object of type \code{clust_irr} 
generated with the function \code{cluster_irr}}
\item{cores}{integer, number of computer cores to use (default = 1)}
\item{as_visnet}{logical, if \code{as_visnet}=TRUE we plot an interactive 
graph with visNetwork. If \code{as_visnet}=FALSE, we plot a static graph 
with igraph.}}
\value{
The output is an igraph plot.

Vertices are clones and edges represent local or global similarities.

The size of the vertices increases linearly as the logarithm of the degree of 
the clonal expansion (number of cells per clone) in the corresponding clones.
}

\description{
This this function creates a joint graph from two or more \code{clust_irr} 
objects, and visualizes the graph.
}

\examples{
# load package input data
data("CDR3ab")
s_1 <- base::data.frame(CDR3b = CDR3ab[1:100, "CDR3b"])
s_2 <- base::data.frame(CDR3b = CDR3ab[101:200, "CDR3b"])
r <- base::data.frame(CDR3b = CDR3ab[1:500, "CDR3b"])


# run 1st analysis -> clust_irr object
o_1 <- cluster_irr(s = s_1, r = r, version = 2, ks = 4)

# run 2nd analysis -> clust_irr object
o_2 <- cluster_irr(s = s_2, r = r, version = 2, ks = 4)

# plot graph with vertices as clones
plot_joint_graph(c(o_1, o_2))
}
