\name{decode_all_communities}
\alias{decode_all_communities}
\title{Decode all graph communities}
\description{
Given a graph based on which we have detected communities (with the function
\code{detect_community}), this function applies \code{decode_communities}
to \bold{all} communities in the graph.  

Each community is partitioned according to user-defined filters on edges 
(\code{edge_filter}) and nodes (\code{node_filter}). This allows extraction 
of subgraphs, cliques, and connected components across the full set of 
communities in the graph.
}
\details{
Internally, this function iterates over all unique community IDs found in the 
vertex attribute \code{community}, and applies \code{decode_communities}.  

\itemize{
\item The \code{edge_filter} controls which edges are retained based on 
their attributes.
\item The \code{node_filter} groups nodes with shared attributes into 
subcomponents.
}
}
\usage{
decode_all_communities(graph, edge_filter, node_filter)
}
\arguments{
\item{graph}{\code{igraph} object that has been analyzed by graph-based 
community detection methods as implemented in \code{detect_community}. 
Must have a vertex attribute named \code{community}.}
\item{edge_filter}{data.frame with edge filters. The data.frame has three 
columns:}
\itemize{
\item name: edge attribute name
\item value: edge attribute value (threshold)
\item operation: logical operation that tells ClustIRR which edge attribute 
values should pass the filter. Possible operations: "<", ">", ">=", "<=", 
"==" and "!=".
}
\item{node_filter}{data.frame with node filters. Groups of nodes that have 
the same attribute values among \bold{ALL} provided attributes will be 
treated as a subcomponent.}
}
\value{
A named list of results, one per community. Each element of the list is the 
output of \code{decode_communities}, containing:
\itemize{
\item community_graph: "filtered" igraph object for the community
\item component_stats: data.frame with summary about each connected component
\item node_summary: data.frame with summary about each node
}
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

# Construct edge and node filters
edge_filter <- rbind(data.frame(name = "nweight", value = 8, operation = ">="),
                     data.frame(name = "ncweight", value = 8, operation = ">="))
node_filter <- data.frame(name = "Ag_gene")

# Decode all communities at once
all_decoded <- decode_all_communities(graph = gcd$graph, 
                                      edge_filter = edge_filter,
                                      node_filter = node_filter)

# Inspect one decoded community
names(all_decoded)
str(all_decoded[[1]]$component_stats)
}
