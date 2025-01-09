\name{detect_communities}
\alias{detect_communities}
\title{Graph-based community detection (GCD)}
\description{
Performs graph-based community detection to find densely connected groups of 
nodes in graph constructed by \code{get_graph} or \code{get_joint_graph}.
}
\details{
ClustIRR employs graph-based community detection (GCD) algorithms, such as 
Louvain or Leiden, to identify densely connected nodes. But first, we must 
decide how to compute a similarity between two nodes, \eqn{i} and \eqn{j}, 
(e.g. TCR clones) based on the similarity scores between their CDR3 sequences 
(computed in \code{clust_irr}) and use this metric as edge weight 
\eqn{\omega(i,j)}.

\bold{Scenario 1} 

If our IRR data data contains CDR3 sequences from only one chain, such as 
CDR3\eqn{\beta}, then \eqn{\omega(i,j)} is defined as

\eqn{\omega(i,j)={\bar{\omega}}^\beta \qquad\text{or}\qquad 
\omega(i,j)={\bar{\omega}}^\beta_c}

The user can decide among the two definitions by specifying 
\itemize{
\item \code{weight} = "ncweight" \eqn{\rightarrow} 
\eqn{\omega(i,j)=\bar{\omega}_c} (default)
\item \code{weight} = "nweight" \eqn{\rightarrow} 
\eqn{\omega(i,j)=\bar{\omega}}
}

\bold{Scenario 2} 

If our IRR data contains CDR3 sequences from both chains (paired data)
To compute the similarity score between TCR clones, \eqn{i} and \eqn{j}, 
we compute the average alignment score (\code{metric}=\bold{average}) 
from their CDR3\eqn{\alpha} and CDR3\eqn{\beta} alignment scores (in the 
next, I will use TCR\eqn{\alpha\beta} as an example, however, this approach 
can also be used to compare TCR\eqn{\gamma\delta} or BCR\emph{IgH-IgL} 
clones): 

\eqn{\omega(i,j)=\dfrac{{\bar{\omega}}^\alpha + {\bar{\omega}}^\beta}{2}
\qquad\text{or}\qquad
\omega(i,j)=\dfrac{{\bar{\omega}}^\alpha_c + {\bar{\omega}}^\beta_c}{2}},

where \eqn{\bar{\omega}^\alpha} and \eqn{\bar{\omega}^\beta} are the
alignment scores for the CDR3\eqn{\alpha} and CDR3\eqn{\beta} sequences,
respectively; and \eqn{\bar{\omega}^\alpha_c} and \eqn{\bar{\omega}^\beta_c} 
are the alignment scores for the CDR3\eqn{\alpha} and CDR3\eqn{\beta} cores,
respectively. Based on this metric, CDR3\eqn{\alpha} and CDR3\eqn{\beta} 
contribute towards the overall similarity of the TCR clones with equal weights.

ClustIRR provides two additional metrics for computing similarity scores 
between TCR clones, including a \code{metric}=\bold{strict}, which assigns 
high similarity score to a pair of TCR clones only if both of their 
CDR3\eqn{\alpha} and CDR3\eqn{\beta} sequence pairs are similar

\eqn{\omega(i,j) = \min({\bar{\omega}}^\alpha, {\bar{\omega}}^\beta)
    \qquad\text{or}\qquad
    \omega(i,j) = \min({\bar{\omega}}^\alpha_c, {\bar{\omega}}^\beta_c)},

and a \code{metric}=\bold{loose}, which assigns high similarity score to
a pair of TCR clones if either of their CDR3\eqn{\alpha} and CDR3\eqn{\beta}
sequence pairs are similar 

\eqn{\omega(i,j) = \max({\bar{\omega}}^\alpha, {\bar{\omega}}^\beta)
    \qquad\text{or}\qquad
    \omega(i,j) = \max({\bar{\omega}}^\alpha_c, {\bar{\omega}}^\beta_c)},
}
\usage{
detect_communities(graph,
                   algorithm = "leiden", 
                   resolution = 1,
                   weight = "ncweight",
                   metric = "average",
                   chains)
}
\arguments{
\item{graph}{\code{igraph} object}
\item{algorithm}{graph-based community detection (GCD) method: leiden (default) 
or louvain.}
\item{resolution}{clustering resolution (default = 1) for the GCD.}
\item{weight}{which edge weight metric (default = ncweight) should be 
used for GCD}
\item{metric}{possible metrics: "average" (default), "strict" or "loose".}
\item{chains}{which chains should be used for clustering? For instance:
chains = "CDR3a"; or chains = CDR3b; or chains = c("CDR3a", "CDR3b").}
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
a <- data.frame(CDR3a = CDR3ab[1:300, "CDR3a"],
                  CDR3b = CDR3ab[1:300, "CDR3b"],
                  clone_size = 1,
                  sample = "a")

b <- data.frame(CDR3a = CDR3ab[201:400, "CDR3a"],
                  CDR3b = CDR3ab[201:400, "CDR3b"],
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
                          weight = "ncweight",
                          metric = "average",
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
