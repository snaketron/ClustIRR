\name{get_community_purity}
\alias{get_community_purity}
\title{Compute community purity with respect to a node feature}
\description{
Given an \code{igraph} object with detected communities (via a node attribute 
called \code{community}) and a node feature (e.g., the output graph of the 
function \code{detect_communities}), this function computes the purity of each 
community with respect to that feature.

If the node feature is numeric (e.g., expression values, clone size, etc.), 
purity is measured by the coefficient of variation (CV). If the node feature 
is categorical (e.g., tissue type, sample name, cell type), purity is measured 
by Gini Impurity (GI).
}
\details{
The function checks for the existence of a \code{community} vertex attribute 
in the graph and for a valid \code{node_feature}.  
Depending on whether the feature is numeric or categorical, it computes CV or
GI, respectively. 

\itemize{
\item For categorical features: communities with low GI contain 
mostly homogeneous feature values.
\item For numeric features: communities with low CV 
are more homogeneous in feature values.
}
}
\usage{
get_community_purity(graph, node_feature)
}
\arguments{
\item{graph}{\code{igraph} object with community membership stored as a 
vertex attribute named \code{community}.}
\item{node_feature}{Character scalar specifying the vertex attribute to 
evaluate (must be either categorical or numeric).}
}
\value{
A \code{data.frame} summarizing purity for each community:
\itemize{
\item For categorical features: columns \code{GI} (Gini impurity), 
\code{n} (community size), \code{community}.
\item For numeric features: columns \code{mean}, \code{sd}, \code{cv} 
(coefficient of variation), \code{n} (community size), \code{community}.
}
}
\examples{
library(igraph)

# create a simple graph
g <- make_ring(10)
V(g)$community <- rep(1:2, each = 5)
V(g)$tissue <- rep(c("A", "B"), times = 5)
V(g)$expression <- rnorm(10, mean = 5, sd = 1)

# categorical feature
get_community_purity(graph = g, node_feature = "tissue")

# numeric feature
get_community_purity(graph = g, node_feature = "expression")
}
