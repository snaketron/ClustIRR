\name{decode_communities}
\alias{decode_communities}
\title{Decode graph communities}
\description{
Given a graph and community ID, the function explores the heterogeneity 
within communities and summarizes its features.
}
\details{

}
\usage{
decode_communities(community_id,
                   graph,
                   edge_attr,
                   node_attr)
}
\arguments{
\item{graph}{\code{igraph} object that has been analyzed by graph-based 
community detection methods as implemented in \code{detect_communities}}
\item{community_id}{which community should be decoded?}
\item{edge_attr}{which edge attributes should be used to decode the community}
\item{iterations}{which node attributes should be used to decode the community}
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

}
