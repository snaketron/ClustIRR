\name{detect_communities}
\alias{detect_communities}
\title{Graph-based community detection (GCD)}
\description{
Graph-based community detection in graphs constructed by \code{get_graph} 
or \code{get_joint_graph}.
}
\details{
ClustIRR employs graph-based community detection (GCD) algorithms, such as 
Louvain, Leiden or InfoMap, to identify \code{communities} of nodes that have 
high density of edges among each other, and low density of edges with nodes 
outside the community.
}
\usage{
detect_communities(graph,
                   weight = "nweight",
                   algorithm = "leiden", 
                   resolution = 1,
                   iterations = 100,
                   chains)
}
\arguments{
\item{graph}{\code{igraph} object}
\item{algorithm}{graph-based community detection (GCD) method: leiden 
(default), louvain or infomap.}
\item{resolution}{clustering resolution (default = 1) for GCD.}
\item{iterations}{clustering iterations (default = 100) for GCD.}
\item{weight}{which edge weight attribute (default = nweight) should be 
used for GCD}
\item{chains}{which chains should be used for clustering? For instance:
chains = "CDR3a"; or chains = "CDR3b"; or chains = c("CDR3a", "CDR3b").}
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
\item{graph}{igraph object, processed graph object.}
\item{graph_structure_quality}{graph modularity and quality (only for Leiden) 
measure of the strength of division of the graph into communities.}
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
c <- clustirr(s = rbind(a, b))

# detect communities
gcd <- detect_communities(graph = c$graph, 
                          algorithm = "leiden",
                          resolution = 1,
                          weight = "ncweight",
                          iterations = 100,
                          chains = c("CDR3a", "CDR3b"))

# look at outputs
names(gcd)

# look at the community occupancymatrix
head(gcd$community_occupancy_matrix)

# look at the community summary
head(gcd$community_summary$wide)

# look at the node summary
head(gcd$node_summary)
}
