\name{detect_communities}
\alias{detect_communities}
\title{Graph-based community detection}
\description{
This algorithm detects densely connected groups of nodes in a immune receptor
repertoire (IRR) specificity graph
}
\usage{
detect_communities(graph,
                   algorithm = "leiden", 
                   resolution = 1,
                   weight_type = "ncweight",
                   chains)
}
\arguments{
\item{graph}{igraph object returned by \code{get_graph} or 
\code{get_joint_graph}}.
\item{algorithm}{graph-based community detection (GCD) algorithm. Two GCD 
algorithms are available: leiden (default) or louvain.}
\item{resolution}{clustering resolution (default = 1) for the GCD algorithm.}
\item{weight_type}{which edge weight metric (default = normalized core weight = 
ncweight) should be used for GCD}
\item{chains}{vector, which chains should be used for clustering? For instance,
chains = "CDR3a", or chains = CDR3b or chains = c("CDR3a", "CDR3b").}
}
\value{
The output is a list with the folling elements:
\item{community_occupancy_matrix}{matrix, rows are communities, columns 
are repertoires, matrix entries are numbers of cells in each community and 
repertoire.}
\item{community_summary}{data.frame, rows are communities and their properties
are provided as columns.}
\item{node_summary}{data.frame, rows are nodes (clones) and their properties 
are provided as columnscontains all user provided.}
\item{graph}{igraph object, processed graph object}
\item{input_config}{list, inputs provided as list.}
}
\examples{
# load package input data
data("CDR3ab", package = "ClustIRR")
a <- data.frame(CDR3a = CDR3ab[1:500, "CDR3a"],
                  CDR3b = CDR3ab[1:500, "CDR3b"],
                  clone_size = 1,
                  sample = "a")

b <- data.frame(CDR3a = CDR3ab[401:900, "CDR3a"],
                  CDR3b = CDR3ab[401:900, "CDR3b"],
                  clone_size = 1,
                  sample = "b")
b$clone_size[1] <- 20

# run ClustIRR analysis
c <- c(cluster_irr(s = a), cluster_irr(s = b))

# get joint graph
jg <- get_joint_graph(clust_irrs = c)

# detect communities
gcd <- detect_communities(graph = jg$graph, 
                          algorithm = "leiden",
                          resolution = 1,
                          weight_type = "ncweight",
                          chains = c("CDR3a", "CDR3b"))

# look at outputs
names(gcd)

# look at the community occupancymatrix
head(gcd$community_occupancy_matrix)

# look at the community summary
head(gcd$community_summary)

# look at the node summary
head(gcd$node_summary)
}
