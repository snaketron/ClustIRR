\name{detect_communities}
\alias{detect_communities}
\title{Graph-based community detection}
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
\item{community_matrix}{matrix, rows are communities, columns are repertoires,
matrix entries are numbers of cells in each community and repertoire.}
\item{community_summary}{data.frame, rows are communities and their properties
are provided as columns.}
\item{node_summary}{data.frame, rows are nodes (clones) and their properties 
are provided as columnscontains all user provided.}
\item{graph}{igraph object, processed graph object}
\item{input_config}{list, inputs provided as list.}
}

\description{
This algorithm finds densely connected groups of nodes in IRR specificity 
graph.
}

\examples{
# load package input data
data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], clone_size = 1)
r <- data.frame(CDR3b = CDR3ab[1:5000, "CDR3b"], clone_size = 1)

# artificially enrich motif 'RWGW' inside sample dataset
substr(x = s$CDR3b[1:20], start = 6, stop = 9) <- "RWGW"

# add an artificial clonal expansion of two sequences to the sample dataset
s <- rbind(s, data.frame(CDR3b = c("CATSRADKPDGLDALETQYF", 
                                   "CATSRAAKPDGLAALSTQYF"),
                         clone_size = 5))

# run ClustIRR analysis
out <- cluster_irr(s = s,
                   r = r,
                   ks = 4,
                   cores = 1,
                   control = list(trim_flank_aa = 3))

# get graph
g <- get_graph(out)

# detect communities
coms <- detect_communities(graph = g$graph, chains = "CDR3b")

# look at outputs
names(coms)

# look at the community matrix
head(coms$community_matrix)

# look at the community summary
head(coms$community_summary)

# look at the node summary
head(coms$node_summary)
}
