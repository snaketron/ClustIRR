\name{decode_communities}
\alias{decode_communities}
\title{Decode graph communities}
\description{
Given a graph based on which we have detected communities (with the function
\code{detect_communities}), and a community ID, the function will try to 
partition the community nodes according to user-defined filters: edge and
node filters.

For instance, the user may only be interested in retaining edges with core 
edge weight > 4; or making sure that nodes that have same 'cell_type' (node 
meta datafrom) are grouped together. Or the user might want to treat all nodes 
that have the same V, D and J gene names and HLA types as subgroups, in which
case all edges between nodes that do not share the same sets of attributes 
are dicarded.

Based on these \code{filters}, ClustIRR will reformat the edges in the 
selected community and then find \bold{connected components} in the 
resulting graph.
}
\details{

}
\usage{
decode_communities(community_id, graph, edge_filter, node_filter)
}
\arguments{
\item{graph}{\code{igraph} object that has been analyzed by graph-based 
community detection methods as implemented in \code{detect_communities}}
\item{community_id}{which community should be decoded?}
\item{edge_filter}{data.frame with edge filters. The deta.frame has three 
columns:}
\itemize{
\item name: edge attribute name
\item value: edge attribute value (threshold)
\item operation: logical operation that tells ClustIRR which edge attribute 
values should pass the filter. Possible operations: "<", ">", ">=", "<=", 
"==" and "!=".
}
\item{node_filter}{a vector with node attributes. Groups of nodes that 
have the same attribute values among \bold{ALL} provided attributes will 
be treated as a subcomponent.}
}
\value{
The output is a "filtered" igraph object.
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

# We "decompose" the communities in the gcd object using decode_community 
# based on the attributes of the edges (edge_filter) and nodes (node_filter).
# We can pick from these edge attributes and create filters:
library(igraph)
edge_attr_names(graph = gcd$graph)

# For instance, the following edge-filter will instruct ClustIRR to keep 
# edges with: edge attributes: nweight>=3 \bold{AND} ncweight>=3
edge_filter <- rbind(data.frame(name = "nweight", value = 3, operation = ">="),
                     data.frame(name = "ncweight", value = 3, operation = ">="))

# In addition, we can construct filters based on the following node attributes:
vertex_attr_names(graph = gcd$graph)

# The following node-filter will instruct ClustIRR to retain edges 
# between nodes that have shared node attributed with respect to ALL 
# of the following node attributes: 
node_filter <- data.frame(name = "Ag_gene")

# Lets inspect community with ID = 1.
c1 <- decode_communities(community_id = 1, 
                         graph = gcd$graph,
                         edge_filter = edge_filter,
                         node_filter = node_filter)

# Plot resulting igraph
par(mar = c(0, 0, 0, 0))
plot(c1, vertex.size = 10)

# Now look at node attributes
as_data_frame(x = c1, what = "vertices")[,c("name", "component_id", 
                                            "CDR3b", "CDR3a", "Ag_gene")]
}

